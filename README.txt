This toolbox provides a MATLAB implementation of several orthogonal
dictionary learning methods for sparse coding. The toolbox focuses
primarily on (but is not limited to) the demonstration of orthogonal
dictionary learning scenarios, where data samples and dictionary atoms can
be represented as square patches. To this end, evaluative and illustrative
features are included, whereas performance might not be cutting-edge.

The following two principal orthognal sparse coding models are covered:

Let X = [x_1, x_2, ..., x_L] be the data matrix, U the dictionary, and A
the coefficient matrix.

(1) the constrained K-sparse model (sparse_mode = 'column_k-sparse'):
  min_{U: U'U=I} ||X - UA||²_F,
      s.t. a_i = arg min_{a: ||a||_0 <= K} ||x_i - Ua_i||²_2 for i=1,...,L

(2) the unconstrained Lagrangian model (sparse_mode = 'hard_thresh'):
  min_{U, A: U'U=I} ||X - UA||²_F + lambda*||A||_0

The key function, orthogonal_dictionary_learning.m, handles the iteration
of learning epochs. The dictionary update depends on the selected method
and is realized by a function handle which is a field of the parameter
struct passed to the function. This function handle is set by the function
default_learning_params.m. orthogonal_dictionary_learning.m is designed to
optionally perform additional computations such as to evaluate the cost
function of the sparse coding model, to calculate the similarity of the
current dictionary with a reference dictionary, to visualize the current
dictionary, or to save the current dictionary as an image file (assuming
that dictionary atoms can be represented as patches).

The supported set of methods

'CA' (--> [1])                          model (1), batch
'OSC' (--> [1])                         model (1), online
'GF-OSC_line_search' (--> [2])          model (1), online
'GF-OSC_cooling_learnrate' (--> [2])    model (1), online
'DDTFC' (--> [3],[4])                   model (2), batch
'lambda-OSC'                            model (2), online
'lambda-GF-OSC_cooling_learnrate'       model (2), online
'lambda-GF-OSC_line_search'             model (2), online

Remember that, depending on the model, params.sparsity_param represents
either the maximal number of non-zero entries per column of A (model (1)),
or the hard threshold such that each entry of A with absolute value below
is set to zero (model (2)).

For a quick start have a look at:

  -->   EXAMPLE_dictionary_recovery.m
        (recovers an orthognal dictionary from snythetic sparse data)

  -->   EXAMPLE_natural_image_patches
        (learns an orthogonal dictionary from patches of an image)

If you have trouble to use the pre-compiled mex files:
 - crop_patches_of_img.mexa64
 - update_lambdaOSC_fast.mexa64
 - update_OSC_fast.mexa64

delete them and compile the mex files on your own by executing the 
commands (within MATLAB and the primary toolbox path)
 - mex crop_patches_of_img.cpp
 - mex update_lambdaOSC_fast.cpp
 - mex update_OSC_fast.cpp

Please report bugs or suggestions for improvements to
Henry.Schuetze@uni-luebeck.de

References:

[1] Schütze, H., Barth, E., and Martinetz, T., "Learning Efficient Data
Representations with Orthogonal Sparse Coding", IEEE Transactions on 
Computational Imaging, 2016 (accepted)

[2] Schütze, H., Barth, E., and Martinetz, T., "Learning Orthogonal Sparse
Representations by Using Geodesic Flow Optimization: International Joint
Conference on Neural Networks (IJCNN), 15540:1-8, 2015

[3] Cai, J.-F., Ji, H., Shen, Z., and Ye G.-B., "Data-driven tight frame
construction and image denoising", Applied and Computational Harmonic
Analysis, vol. 37, no. 1, pp. 89--105, 2014

[4] Bao, C., Cai, J.-F., and Ji, H., "Fast sparsity-based orthogonal dictionary
learning for image restoration", IEEE International Conference on Computer
Vision (ICCV), 2013