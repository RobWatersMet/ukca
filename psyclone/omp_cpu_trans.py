from typing import List, Tuple

from psyclone.psyir import nodes
from psyclone.psyir.transformations.arrayassignment2loops_trans import (
    ArrayAssignment2LoopsTrans,
)
# from psyclone.psyir.transformations.parallel_loop_trans import ParallelLoopTrans # TODO
from psyclone.transformations import TransformationError, OMPLoopTrans
# from psytran.clauses import has_collapse_clause
from psytran.convert import convert_array_notation
from psytran.directives import (
     apply_loop_directive,
)
# from psytran.family import has_ancestor
# from psytran.loop import is_parallelisable, is_perfectly_nested


# Example excluded loops
# TODO - need to have if bounds is Idx - i.e. psyclone added it, dont have to specify the number
EXCLUDED_LOOPS = {
    'asad_sparse_vars.F90': {
        'spfuljac': [('jc', 'nstat'), ('j3', 'nmsjac'), ('jn', 'jpmsp'),
                     ('jn', None)]
    }
}


EXCLUDED_NODE_TYPES = (nodes.Call, nodes.CodeBlock, nodes.Return)
EXCLUDED_VARIABLES = ("jit", "iter")


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
    print(f"Applying OpenMP transformations to routine: {routine.name}")

    def is_problematic(node):
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
            [ref.name in EXCLUDED_VARIABLES for ref in
                node.walk(nodes.Reference)]
        )
        return excluded or bad_ref

    def has_excluded_loop(node, excluded_loops):

        if excluded_loops is None:
            return False
        else:
            return _contains_excluded_loop(node, excluded_loops)

    # Normalise all loops within the schedule
    convert_array_notation(routine)
    for assign in routine.walk(nodes.Assignment):
        try:
            ArrayAssignment2LoopsTrans().apply(assign)
        except TransformationError:
            pass


    blocks = routine.get_sibling_lists(nodes.Node)
    by_depth = {}

    def add_block(block):
        """
        Check that a block of nodes is valid and add it as a candidate to apply a
        kernels directive to.
        """
        if not (block and any([isinstance(node, nodes.Loop) for node in block])):
            return
        depth = block[0].depth
        if depth not in by_depth:
            by_depth[depth] = []
        by_depth[depth].append(block)

    # Drop problematic nodes and blocks
    for block in blocks:
        local_list = []
        for node in block:
            if is_problematic(node) or has_excluded_loop(node, excluded_loops):
                add_block(local_list)
                local_list = []
            else:
                local_list.append(node)
        add_block(local_list)
        depth = block[0].depth
        if depth in by_depth and not by_depth[depth]:
            by_depth.pop(depth)


    # TODO: Develop OpenMP transformations, noting kernels not valid

    #     # Extract blocks of potentially parallelisable code
    #     blocks = schedule.get_sibling_lists(nodes.Node)
    #     by_depth = {}

    #     def add_block(block):
    #         """
    #         Check that a block of nodes is valid and add it as a candidate to apply a
    #         kernels directive to.
    #         """
    #         if not (block and any([isinstance(node, nodes.Loop) for node in block])):
    #             return
    #         depth = block[0].depth
    #         if depth not in by_depth:
    #             by_depth[depth] = []
    #         by_depth[depth].append(block)

    #     # Drop problematic nodes and blocks
    #     for block in blocks:
    #         local_list = []
    #         for node in block:
    #             if is_problematic(node):
    #                 add_block(local_list)
    #                 local_list = []
    #             else:
    #                 local_list.append(node)
    #         add_block(local_list)
    #         depth = block[0].depth
    #         if depth in by_depth and not by_depth[depth]:
    #             by_depth.pop(depth)

    #     # Apply OpenACC kernels directives
    #     for blocks in by_depth.values():
    #         for block in blocks:
    #             if not has_kernels_directive(block):
    #                 local_list = []
    #                 if any([isinstance(node, nodes.IfBlock) for node in block]):
    #                     for node in block:
    #                         if isinstance(node, nodes.IfBlock):
    #                             apply_kernels_dir(local_list)
    #                             apply_loop_trans_block(local_list)
    #                             local_list = []
    #                         else:
    #                             local_list.append(node)
    #                     block = local_list
    #                 apply_kernels_dir(block)
    #             apply_loop_trans_block(block)

    return


def _contains_excluded_loop(node: nodes.Node, excluded_loops: List[Tuple]
                            ) -> bool:
    """
    Private function to work out if a node has an excluded loop as a child
    within it.

    Args:
        node (nodes.Node): the pysir node to look below
        excluded_loops (List[Tuple]): list of excluded loops in format -
            (loop index variable, upper bounds variable)

    Returns:
        bool: whether and excluded loop is a child of the nde
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
            if excluded_loop[0] == variable and \
                    excluded_loop[1] in references:
                return True

    return False


# def apply_kernels_dir(block):
#     """
#     Apply OpenACC kernels directives to a block of code.
#     """
#     if not block:
#         return
#     try:
#         apply_kernels_directive(block, options={"force": True})
#     except TransformationError:
#         pass
#
#
# def apply_loop_trans(loop):
#     """
#     Apply OpenACC loop directives to the loop.
#
#     :arg loop: the :class:`nodes.Loop` instance to apply `!$acc loop` to
#     """
#     assert isinstance(loop, nodes.Loop)
#     if not has_kernels_directive(loop):
#         apply_kernels_directive(loop)
#     if has_collapse_clause(loop):
#         return
#     parallelisable = is_parallelisable(loop)
#     perfect_nest = (
#         is_perfectly_nested(loop)
#         and not has_ancestor(loop, name="jl")
#     )
#
#     # Determine clauses to apply
#     options = {"force": True}
#     if perfect_nest and parallelisable:
#         options["gang"] = True
#         options["vector"] = True
#         options["collapse"] = True
#     elif loop.variable.name == "jl" and parallelisable:
#         options["gang"] = True
#         options["vector"] = True
#     else:
#         options["sequential"] = True
#
#     apply_loop_directive(loop, options=options)
#
#
# def apply_loop_trans_block(block):
#     """
#     Apply OpenACC loop directives to all loops within a block of code.
#     """
#     for node in block:
#         if isinstance(node, nodes.Loop) and not has_loop_directive(node):
#             apply_loop_trans(node)
