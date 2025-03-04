library(scistreer)

test_that("Small input works", {
    tree = run_scistree(P_small, verbose = F)
    expect_equal(round(tree$likelihood, 3), -1606.497)
})

test_that("Score tree works", {
    l = score_tree(tree_small, P_small)$l_tree
    expect_equal(round(l, 3), -1606.497)
})

test_that("Score tree C++ works", {
    tree_small = reorder(tree_small, order = 'postorder')
    l = score_tree_cpp(tree_small$edge, P_small)
    expect_equal(round(l, 3), -1606.497)
})

# test_that("Conversion between phylo and tree graph works", {
#     tree_small_new = to_phylo(annotate_tree(tree_small, P_small))
#     l0 = score_tree(to_phylo(gtree_small), P_small)$l_tree
#     l1 = score_tree(tree_small, P_small)$l_tree
# })