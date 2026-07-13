// Algebra oracles mirroring main/joint_meta.cpp inv_update_* (incl. denom / cond guards).

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>

namespace {

using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
using SpMat = Eigen::SparseMatrix<double>;

bool inv_update_forward_append(const Mat& R_inv, const Vec& c, double d, Mat& R_inv_out) {
    const int k = static_cast<int>(R_inv.rows());
    if (k == 0) {
        if (!(d > 0.0)) return false;
        R_inv_out.resize(1, 1);
        R_inv_out(0, 0) = 1.0 / d;
        return true;
    }
    Vec v = R_inv * c;
    const double denom = d - c.dot(v);
    if (!(denom > 0.0)) return false;
    const double s = 1.0 / denom;
    R_inv_out.resize(k + 1, k + 1);
    R_inv_out.topLeftCorner(k, k) = R_inv + s * v * v.transpose();
    R_inv_out.topRightCorner(k, 1) = -s * v;
    R_inv_out.bottomLeftCorner(1, k) = -s * v.transpose();
    R_inv_out(k, k) = s;
    return true;
}

Mat update_inv_backward(const Mat& R_inv_pre, int remove_index) {
    Mat R_inv_post =
        R_inv_pre - R_inv_pre.col(remove_index) * R_inv_pre.row(remove_index) /
                        R_inv_pre(remove_index, remove_index);
    const int dim = static_cast<int>(R_inv_post.rows());
    R_inv_post.block(remove_index, 0, dim - remove_index - 1, dim) =
        R_inv_post.block(remove_index + 1, 0, dim - remove_index - 1, dim);
    R_inv_post.conservativeResize(dim - 1, Eigen::NoChange);
    R_inv_post.block(0, remove_index, dim - 1, dim - remove_index - 1) =
        R_inv_post.block(0, remove_index + 1, dim - 1, dim - remove_index - 1);
    R_inv_post.conservativeResize(Eigen::NoChange, dim - 1);
    return R_inv_post;
}

Mat equicorr3(double rho) {
    Mat R(3, 3);
    R.setConstant(rho);
    R.diagonal().setOnes();
    return R;
}

double max_collinear(const Mat& R_inv, const Vec& diagB) {
    return (1.0 - Vec::Ones(R_inv.rows()).array() / (diagB.array() * R_inv.diagonal().array()))
        .maxCoeff();
}

}  // namespace

TEST(CojoAlgebra, ForwardInverseUpdateMatchesDirectInverse) {
    Mat R_pre(2, 2);
    R_pre << 1.0, 0.2, 0.2, 1.0;

    Vec c(2);
    c << 0.3, -0.1;
    const double d = 1.0;

    Mat R_post(3, 3);
    R_post.topLeftCorner(2, 2) = R_pre;
    R_post.topRightCorner(2, 1) = c;
    R_post.bottomLeftCorner(1, 2) = c.transpose();
    R_post(2, 2) = d;

    Mat R_inv_pre = R_pre.ldlt().solve(Mat::Identity(2, 2));
    Mat R_inv_expected = R_post.ldlt().solve(Mat::Identity(3, 3));
    Mat R_inv_post;
    ASSERT_TRUE(inv_update_forward_append(R_inv_pre, c, d, R_inv_post));

    EXPECT_TRUE(R_inv_post.isApprox(R_inv_expected, 1e-12));
}

TEST(CojoAlgebra, ForwardAppendRejectsNonPositiveSchur) {
    Mat R_pre(2, 2);
    R_pre << 1.0, 0.9, 0.9, 1.0;
    Mat R_inv_pre = R_pre.ldlt().solve(Mat::Identity(2, 2));

    // c identical to an existing column → Schur complement ≈ 0 / negative.
    Vec c(2);
    c << 1.0, 0.9;
    Mat R_inv_out;
    EXPECT_FALSE(inv_update_forward_append(R_inv_pre, c, 1.0, R_inv_out));
}

TEST(CojoAlgebra, BackwardInverseUpdateMatchesDirectInverse) {
    Mat R(3, 3);
    R << 1.0, 0.2, -0.1, 0.2, 1.0, 0.3, -0.1, 0.3, 1.0;

    Mat R_inv = R.ldlt().solve(Mat::Identity(3, 3));

    const int remove_index = 1;
    Mat R_reduced(2, 2);
    R_reduced << 1.0, -0.1, -0.1, 1.0;

    Mat R_inv_expected = R_reduced.ldlt().solve(Mat::Identity(2, 2));
    Mat R_inv_post = update_inv_backward(R_inv, remove_index);

    EXPECT_TRUE(R_inv_post.isApprox(R_inv_expected, 1e-12));
}

TEST(CojoAlgebra, CondGuardEstimateMatchesStockThreshold) {
    // Stock uses sqrt(max(D)/min(D)) > 30 on SimplicialLDLT of _B.
    // For equicorrelation this only exceeds 30 at extreme ρ, where the
    // per-SNP collinearity measure also exceeds 0.9 — so the two guards
    // largely overlap on correlation-scale blocks. Still assert the D-ratio path.
    const double rho = 0.9999;
    Mat R = equicorr3(rho);
    SpMat Rs = R.sparseView();
    Eigen::SimplicialLDLT<SpMat> ldlt(Rs);
    ASSERT_EQ(ldlt.info(), Eigen::Success);

    const double dmin = ldlt.vectorD().minCoeff();
    const double dmax = ldlt.vectorD().maxCoeff();
    ASSERT_GT(dmin, 0.0);
    const double cond_est = std::sqrt(dmax / dmin);
    EXPECT_GT(cond_est, 30.0);

    Mat R_inv = ldlt.solve(Mat::Identity(3, 3));
    Vec diagB = Vec::Ones(3);
    EXPECT_GT(max_collinear(R_inv, diagB), 0.9);
}

TEST(CojoAlgebra, NearDuplicateRejectedByGuards) {
    Mat R2(2, 2);
    R2 << 1.0, 0.4, 0.4, 1.0;
    Mat R2_inv = R2.ldlt().solve(Mat::Identity(2, 2));

    Vec c(2);
    c << 0.999, 0.4;  // nearly copy of SNP0
    const double d = 1.0;

    Mat R3(3, 3);
    R3.topLeftCorner(2, 2) = R2;
    R3.topRightCorner(2, 1) = c;
    R3.bottomLeftCorner(1, 2) = c.transpose();
    R3(2, 2) = d;

    Mat R3_inv_incr;
    const bool schur_ok = inv_update_forward_append(R2_inv, c, d, R3_inv_incr);

    SpMat R3s = R3.sparseView();
    Eigen::SimplicialLDLT<SpMat> ldlt(R3s);
    const bool ldlt_ok = ldlt.info() == Eigen::Success && ldlt.vectorD().minCoeff() > 0;
    bool cond_ok = false;
    bool col_ok = false;
    if (ldlt_ok) {
        const double est =
            std::sqrt(ldlt.vectorD().maxCoeff() / ldlt.vectorD().minCoeff());
        cond_ok = est <= 30.0;
        Mat R3_inv = ldlt.solve(Mat::Identity(3, 3));
        col_ok = max_collinear(R3_inv, Vec::Ones(3)) <= 0.9;
    }
    const bool accept_full = ldlt_ok && cond_ok && col_ok;
    EXPECT_FALSE(accept_full);
    // Schur alone may or may not reject; full guard set must.
    (void)schur_ok;
}

// Moving last Z-cache row to sorted insert position (append-then-permute).
TEST(CojoAlgebra, MoveLastRowToPosPreservesOtherRows) {
    Mat Z(3, 4);
    Z << 1, 2, 3, 4,
         5, 6, 7, 8,
         9, 10, 11, 12;
    // Simulate append of row [90,91,92,93] then move to pos=1
    Mat Z2(4, 4);
    Z2.topRows(3) = Z;
    Z2.row(3) << 90, 91, 92, 93;
    const int sorted_pos = 1;
    const int k = 3; // last index before permute (= new_rows-1)
    Mat out = Z2;
    Eigen::RowVectorXd z_last = out.row(k);
    for (int r = k; r > sorted_pos; r--) out.row(r) = out.row(r - 1);
    out.row(sorted_pos) = z_last;

    EXPECT_DOUBLE_EQ(out(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(out(1, 0), 90.0);
    EXPECT_DOUBLE_EQ(out(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(out(3, 0), 9.0);
}

// Centered product from dosage sums matches explicit float center (GCTA /n scale).
TEST(CojoAlgebra, BitStyleCenteredProductMatchesFloat) {
    // Two SNPs, 8 individuals, no missing.
    const int n = 8;
    const double g1[8] = {0, 1, 2, 1, 0, 2, 1, 0};
    const double g2[8] = {2, 2, 1, 0, 1, 1, 0, 2};
    double mu1 = 0, mu2 = 0;
    for (int i = 0; i < n; i++) {
        mu1 += g1[i];
        mu2 += g2[i];
    }
    mu1 /= n;
    mu2 /= n;
    double float_dot = 0;
    for (int i = 0; i < n; i++) float_dot += (g1[i] - mu1) * (g2[i] - mu2);
    float_dot /= n;

    int S12 = 0, S1 = 0, S2 = 0, N_both = n;
    for (int i = 0; i < n; i++) {
        S12 += (int)(g1[i] * g2[i]);
        S1 += (int)g1[i];
        S2 += (int)g2[i];
    }
    const double bits_dot =
        (S12 - mu1 * S2 - mu2 * S1 + (double)N_both * mu1 * mu2) / (double)n;
    EXPECT_NEAR(bits_dot, float_dot, 1e-12);
}
