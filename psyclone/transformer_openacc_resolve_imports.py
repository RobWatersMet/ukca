from typing import List, Tuple, Optional
from psyclone.psyir import nodes
from psyclone.psyir.transformations.arrayassignment2loops_trans import (
    ArrayAssignment2LoopsTrans,
)
from psyclone.transformations import TransformationError
from psyacc import (
    convert_array_notation,
    has_kernels_directive,
    apply_kernels_directive,
    has_collapse_clause,
    is_parallelisable,
    is_perfectly_nested,
    has_ancestor,
    apply_loop_directive,
    has_loop_directive
)


RESOLVE_IMPORTS = True
FILES_TO_SKIP = ['ukca_cspecies.F90', 'ukca_top_boundary.F90',
                 'ukca_deriv_aero.F90', 'ukca_deriv.F90',
                 'ukca_chem_offline.F90', 'ukca_chem_raqaero_mod.F90',
                 'get_molmass_mod.F90', 'ukca_chem_diags_mod.F90',
                 'ukca_chem1_dat.F90', 'ukca_chem_master.F90',
                 'ukca_chem_raq.F90', 'ukca_chem_aer.F90',
                 'ukca_aerod.F90', 'ukca_surfddr.F90']

EXCLUDED_VARIABLES = ("jit", "iter", "advt")
EXCLUDED_NODE_TYPES = (nodes.Call, nodes.CodeBlock, nodes.Return)


EXCLUDED_LOOPS = {
    'asad_sparse_vars.F90': {
        'spfuljac': [('jc', 'nstat'), ('j3', 'nmsjac'), ('jn', 'jpmsp'),
                     ('jn', None)]
    }
}

# TODO - do not use kernel regions for assignments


def trans(psy):

    print("Applying OpenACC kernels transformations...")

    for routine in psy.walk(nodes.Routine):

        excluded_loops = None
        if psy.name in EXCLUDED_LOOPS.keys():
            if routine.name in EXCLUDED_LOOPS[psy.name].keys():
                excluded_loops = EXCLUDED_LOOPS[psy.name][routine.name]

        tranform_routine(routine, excluded_loops)

    return psy


def tranform_routine(routine: nodes.Routine, excluded_loops: List[Tuple]):

    # Normalise all loops within the schedule
    convert_array_notation(routine)
    for assign in routine.walk(nodes.Assignment):

        try:
            ArrayAssignment2LoopsTrans().apply(assign)
        except TransformationError:
            pass

    # Extract blocks of potentially parallelisable code
    blocks = routine.get_sibling_lists(nodes.Node)
    by_depth = {}

    def add_block(block):
        """
        Check that a block of nodes is valid and add it as a candidate to
        apply a kernels directive to.
        """
        if not (block and any([isinstance(node, nodes.Loop) for node in
                               block])):
            return
        depth = block[0].depth
        if depth not in by_depth:
            by_depth[depth] = []
        by_depth[depth].append(block)

    # Drop problematic nodes and blocks
    for block in blocks:
        local_list = []
        for node in block:
            if is_problematic(node, excluded_loops):
                add_block(local_list)
                local_list = []
            else:
                local_list.append(node)
        add_block(local_list)
        depth = block[0].depth
        if depth in by_depth and not by_depth[depth]:
            by_depth.pop(depth)

    # Apply OpenACC kernels directives
    for blocks in by_depth.values():
        for block in blocks:
            if not has_kernels_directive(block):
                local_list = []
                if any([isinstance(node, nodes.IfBlock) for node in block]):
                    for node in block:
                        if isinstance(node, nodes.IfBlock):
                            apply_kernels_dir(local_list)
                            apply_loop_trans_block(local_list)
                            local_list = []
                        else:
                            local_list.append(node)
                    block = local_list
                apply_kernels_dir(block)
            apply_loop_trans_block(block)

    return


def is_problematic(node: nodes.Node,
                   excluded_loops: Optional[List[Tuple]] = None):
    """
    Determine whether a Node contains excluded node types or variables.
    """
    intrinsic_calls = node.walk(nodes.IntrinsicCall)
    excluded = any(
        [
            call
            for call in node.walk(EXCLUDED_NODE_TYPES)
            if call not in intrinsic_calls
        ]
    )
    bad_ref = any(
        [ref.name in EXCLUDED_VARIABLES for ref in node.walk(
            nodes.Reference)]
    )

    if excluded_loops is None:
        has_excluded_loop = False
    else:
        has_excluded_loop = contains_excluded_loop(node, excluded_loops)

    return excluded or bad_ref or has_excluded_loop


def contains_excluded_loop(node: nodes.Node, excluded_loops: List[Tuple]
                           ) -> bool:
    """
    _summary_

    Args:
        node (nodes.Node): _description_
        excluded_loops (List[Tuple]): _description_

    Returns:
        bool: _description_
    """
    for loop in node.walk(nodes.Loop):
        if any([type(i) is nodes.array_reference.ArrayReference for i in
                loop.children]):
            references = [i.name for i in loop.children if type(i) is
                          nodes.array_reference.ArrayReference]
        elif any([type(i) is nodes.reference.Reference for i in
                  loop.children]):
            references = [i.name for i in loop.children if type(i) is
                          nodes.reference.Reference]
        else:
            references = [None]

        variable = loop.variable.name

        for excluded_loop in excluded_loops:
            try:
                if excluded_loop[0] == variable and excluded_loop[1] in references:
                    return True
            except TypeError:
                import pdb
                pdb.set_trace()
    return False


def apply_kernels_dir(block):
    """
    Apply OpenACC kernels directives to a block of code.
    """
    if not block:
        return
    try:
        apply_kernels_directive(block, options={"force": True})
    except TransformationError:
        pass


def apply_loop_trans(loop):
    """
    Apply OpenACC loop directives to the loop.

    :arg loop: the :class:`nodes.Loop` instance to apply `!$acc loop` to
    """
    assert isinstance(loop, nodes.Loop)
    if not has_kernels_directive(loop):
        try:
            apply_kernels_directive(loop)
        except TransformationError as err:
            print('RW - unable to parallelise loop with variable'
                  f' {loop.variable.name} at position {loop.position}')
            print(f'RW - reason: {err}')
            return
    if has_collapse_clause(loop):
        return
    try:
        parallelisable = is_parallelisable(loop)
    except NotImplementedError as err:
        print('RW - unable to parallelise loop with variable '
              f'{loop.variable.name} at position {loop.position}')
        print(f'RW - reason: {err}')
        return

    perfect_nest = (
        is_perfectly_nested(loop)
        and not has_ancestor(loop, name="jl")
    )

    # Determine clauses to apply
    options = {"force": True}
    if perfect_nest and parallelisable:
        options["gang"] = True
        options["vector"] = True
        nloops = len({node.depth for node in loop.walk(nodes.Loop)})
        if nloops > 1:
            options["collapse"] = nloops
    elif loop.variable.name == "jl" and parallelisable:
        options["gang"] = True
        options["vector"] = True
    else:
        options["sequential"] = True

    apply_loop_directive(loop, options=options)


def apply_loop_trans_block(block):
    """
    Apply OpenACC loop directives to all loops within a block of code.
    """
    for node in block:
        if isinstance(node, nodes.Loop) and not has_loop_directive(node):
            apply_loop_trans(node)
