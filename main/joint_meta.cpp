/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for GWAS simulation
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "gcta.h"
#include <algorithm>
#include <vector>

namespace {

// Rank-1 block-inverse update when appending a row/col to an SPD matrix.
// R_new = [ R , c ; c' , d ], given R_inv = R^{-1}.
// Returns false if the Schur complement is non-positive (non-PD / singular).
bool inv_update_forward_append(const eigenMatrix &R_inv, const eigenVector &c, double d,
                               eigenMatrix &R_inv_out)
{
    const int k = static_cast<int>(R_inv.rows());
    if (k == 0) {
        if (!(d > 0.0)) return false;
        R_inv_out.resize(1, 1);
        R_inv_out(0, 0) = 1.0 / d;
        return true;
    }
    eigenVector v = R_inv * c;
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

// Move last row/col to position `pos` (used after append-at-end to match sorted SNP order).
eigenMatrix move_last_to_pos(const eigenMatrix &M, int pos)
{
    const int n = static_cast<int>(M.rows());
    if (pos < 0 || pos >= n || pos == n - 1) return M;
    eigenMatrix R(n, n);
    // new row i comes from old index src[i]
    for (int i = 0; i < n; i++) {
        const int si = (i < pos) ? i : (i == pos ? n - 1 : i - 1);
        for (int j = 0; j < n; j++) {
            const int sj = (j < pos) ? j : (j == pos ? n - 1 : j - 1);
            R(i, j) = M(si, sj);
        }
    }
    return R;
}

// Rank-1 downdate removing row/col `pos`.
void inv_update_backward(const eigenMatrix &R_inv, int pos, eigenMatrix &R_inv_out)
{
    const int n = static_cast<int>(R_inv.rows());
    if (n <= 1) {
        R_inv_out.resize(0, 0);
        return;
    }
    eigenMatrix tmp =
        R_inv - R_inv.col(pos) * R_inv.row(pos) / R_inv(pos, pos);
    R_inv_out.resize(n - 1, n - 1);
    int ii = 0;
    for (int i = 0; i < n; i++) {
        if (i == pos) continue;
        int jj = 0;
        for (int j = 0; j < n; j++) {
            if (j == pos) continue;
            R_inv_out(ii, jj) = tmp(i, j);
            jj++;
        }
        ii++;
    }
}

// Extract off-diagonal column of lower-triangular sparse B at `pos`, skipping the diagonal,
// in order of the other indices (i.e. old sorted order).
eigenVector extract_sparse_col_excl_diag(const eigenSparseMat &B, int pos, int n)
{
    eigenVector c(n - 1);
    int t = 0;
    for (int i = 0; i < n; i++) {
        if (i == pos) continue;
        if (i > pos)
            c[t] = B.coeff(i, pos);
        else
            c[t] = B.coeff(pos, i);
        t++;
    }
    return c;
}

}  // namespace

void gcta::set_diff_freq(double freq_diff){
    _diff_freq = freq_diff;
}
void gcta::read_metafile(string metafile, bool GC, double GC_val) {
    double freq_diff_thresh = _diff_freq;
    LOGGER << "\nReading GWAS summary-level statistics from [" + metafile + "] ..." << endl;
    if (_bim_region_filter) {
        LOGGER << "BIM region filter active (chr=" << _bim_region_chr
               << " bp=[" << _bim_region_start << "," << _bim_region_end
               << "]): matching .ma SNPs to region index only." << endl;
    }
    ifstream Meta(metafile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + metafile + "] to read.");

    int i = 0, count = 0;
    double f_buf = 0.0, b_buf = 0.0, se_buf = 0.0, p_buf = 0.0, N_buf = 0.0, Vp_buf = 0.0, GC_buf = 0.0, chi_buf = 0.0, h_buf = 0.0;
    string A1_buf, A2_buf;
    string snp_buf, str_buf0, str_buf;

    vector<string> snplist, vs_buf, bad_snp;
    vector<string> ref_A1_buf, ref_A2_buf, bad_A1, bad_A2, bad_refA;
    vector<double> freq_buf, beta_buf, beta_se_buf, pval_buf, N_o_buf, Vp_v_buf, GC_v_buf;
    map<string, int>::iterator iter;
    getline(Meta, str_buf); // the header line
    if (StrFunc::split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in the input file [" + metafile + "].");
    _jma_Vp = 0.0;
    _GC_val = -1;
    int matched_ma = 0, skipped_ma = 0;
    while (Meta) {
        getline(Meta, str_buf0);
        stringstream iss(str_buf0);
        iss >> snp_buf >> A1_buf >> A2_buf;
        // With a region-filtered BIM, most .ma lines miss the map — check before allele work.
        iter = _snp_name_map.find(snp_buf);
        const bool in_region_index = (iter != _snp_name_map.end());
        if (_bim_region_filter && !in_region_index) {
            // Still need Vp from all lines for phenotypic variance median (parity with stock).
            iss >> str_buf;
            f_buf = atof(str_buf.c_str());
            iss >> str_buf;
            if (str_buf == "NA" || str_buf == ".") { skipped_ma++; continue; }
            b_buf = atof(str_buf.c_str());
            iss >> str_buf;
            if (str_buf == "NA" || str_buf == "." || str_buf == "0") { skipped_ma++; continue; }
            se_buf = atof(str_buf.c_str());
            iss >> str_buf;
            if (str_buf == "NA" || str_buf == ".") { skipped_ma++; continue; }
            p_buf = atof(str_buf.c_str());
            iss >> str_buf;
            if (str_buf == "NA" || str_buf == ".") { skipped_ma++; continue; }
            N_buf = atof(str_buf.c_str());
            if (Meta.eof()) break;
            if (N_buf < 10) LOGGER.e(0, "invalid sample size in line:\n\"" + str_buf0 + "\"");
            h_buf = 2.0 * f_buf * (1.0 - f_buf);
            Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
            if (Vp_buf < 0.0) LOGGER.e(0, "in line:\n\"" + str_buf0 + "\"");
            Vp_v_buf.push_back(Vp_buf);
            if (GC) {
                GC_buf = b_buf * b_buf / se_buf / se_buf;
                if (GC_buf < 0) LOGGER.e(0, "in line:\n\"" + str_buf0 + "\"");
                GC_v_buf.push_back(GC_buf);
            }
            count++;
            skipped_ma++;
            continue;
        }
        StrFunc::to_upper(A1_buf);
        StrFunc::to_upper(A2_buf);
        iss >> str_buf;
        f_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        b_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == "." || str_buf == "0") continue;
        se_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        p_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        N_buf = atof(str_buf.c_str());
        if (N_buf < 10) LOGGER.e(0, "invalid sample size in line:\n\"" + str_buf0 + "\"");
        if (Meta.eof()) break;
        iter = _snp_name_map.find(snp_buf);
        h_buf = 2.0 * f_buf * (1.0 - f_buf);
        Vp_buf = h_buf * N_buf * se_buf * se_buf + h_buf * b_buf * b_buf * N_buf / (N_buf - 1.0);
        if (Vp_buf < 0.0) LOGGER.e(0, "in line:\n\"" + str_buf0 + "\"");
        Vp_v_buf.push_back(Vp_buf);
        if (GC) {
            GC_buf = b_buf * b_buf / se_buf / se_buf;
            if (GC_buf < 0) LOGGER.e(0, "in line:\n\"" + str_buf0 + "\"");
            GC_v_buf.push_back(GC_buf);
        }
        count++;
        if (iter == _snp_name_map.end()) { skipped_ma++; continue; }
        i = iter->second;
        if (A1_buf != _allele1[i] && A1_buf != _allele2[i]) {
            bad_snp.push_back(_snp_name[i]);
            bad_A1.push_back(_allele1[i]);
            bad_A2.push_back(_allele2[i]);
            bad_refA.push_back(A1_buf);
            continue;
        }
        snplist.push_back(snp_buf);
        ref_A1_buf.push_back(A1_buf);
        ref_A2_buf.push_back(A2_buf);
        freq_buf.push_back(f_buf);
        beta_buf.push_back(b_buf);
        beta_se_buf.push_back(se_buf);
        pval_buf.push_back(p_buf);
        N_o_buf.push_back(N_buf);
        matched_ma++;
    }
    Meta.close();
    LOGGER << "GWAS summary statistics of " << count << " SNPs read from [" + metafile + "]." << endl;
    if (_bim_region_filter) {
        LOGGER << "MA region match: " << matched_ma << " in index, " << skipped_ma << " outside/unmatched." << endl;
    }
    _jma_Vp = CommFunc::median(Vp_v_buf);
    LOGGER << "Phenotypic variance estimated from summary statistics of all " << count << " SNPs: " << _jma_Vp << " (variance of logit for case-control studies)." << endl;
    if (GC) {
        if (GC_val > 0) {
            _GC_val = GC_val;
            LOGGER << "User specified genomic inflation factor: " << _GC_val << endl;
        } else {
            _GC_val = CommFunc::median(GC_v_buf) / 0.455;
            LOGGER << "Genomic inflation factor calculated from " << count << " SNPs: " << _GC_val << endl;
        }
        LOGGER << "p-values will be adjusted by the genomic control approach." << endl;
    }

    LOGGER << "Matching the GWAS meta-analysis results to the genotype data ..." << endl;
    update_id_map_kp(snplist, _snp_name_map, _include);
    if (_mu.empty()) calcu_mu();

    vector<int> indx;
    vector<string> snplist_freq;
    vector<string> bad_snp_freq, bad_A1_freq, bad_A2_freq, bad_refA_freq, bad_otherA_freq;
    vector<double> bad_freq_value, bad_freq_ma;
    map<string, int> id_map;
    for (i = 0; i < snplist.size(); i++) id_map.insert(pair<string, int>(snplist[i], i));
    for (i = 0; i < _include.size(); i++) {
        int include_i = _include[i];
        bool flip_flag = false;
        string cur_snp_name = _snp_name[include_i];
        iter = id_map.find(cur_snp_name);
        _ref_A[include_i] = ref_A1_buf[iter->second];
        _other_A[include_i] = ref_A2_buf[iter->second];
        if (!_mu.empty() && ref_A1_buf[iter->second] == _allele2[include_i]) {
            _mu[include_i] = 2.0 - _mu[include_i];
            flip_flag = true;
        }
        double cur_freq_value =  _mu[include_i] / 2.0;
        double freq_diff = abs(cur_freq_value - freq_buf[iter->second]);
        if(freq_diff < freq_diff_thresh){
            snplist_freq.push_back(cur_snp_name);
            indx.push_back(iter->second);
        }else{
            if(flip_flag) cur_freq_value = 1.0 - cur_freq_value;
            bad_snp_freq.push_back(cur_snp_name);
            bad_A1_freq.push_back(_allele1[include_i]);
            bad_A2_freq.push_back(_allele2[include_i]);
            bad_refA_freq.push_back(_ref_A[include_i]);
            bad_otherA_freq.push_back(_other_A[include_i]);
            bad_freq_value.push_back(cur_freq_value);
            bad_freq_ma.push_back(freq_buf[iter->second]);
        }
        
    }
    update_id_map_kp(snplist_freq, _snp_name_map, _include);

    if (!bad_snp.empty()) {
        string badsnpfile = _out + ".badsnps";
        ofstream obadsnp(badsnpfile.c_str());
        obadsnp << "SNP\tA1\tA2\tRefA" << endl;
        for (i = 0; i < bad_snp.size(); i++) obadsnp << bad_snp[i] << "\t" << bad_A1[i] << "\t" << bad_A2[i] << "\t" << bad_refA[i] << endl;
        obadsnp.close();
        LOGGER << "Warning: can't match the reference alleles of " << bad_snp.size() << " SNPs to those in the genotype data. These SNPs have been saved in [" + badsnpfile + "]." << endl;
    }
    if(!bad_snp_freq.empty()){
        string badsnpfile = _out + ".freq.badsnps";
        ofstream obadsnp(badsnpfile.c_str());
        obadsnp << "SNP\tsummary_A1\tsummary_A2\tsummary_A1_freq\tgeno_A1\tgeno_A2\tgeno_A1_freq" << endl;
        for (int i = 0; i < bad_snp_freq.size(); i++){
            obadsnp << bad_snp_freq[i] << "\t"
                    << bad_refA_freq[i] << "\t"
                    << bad_otherA_freq[i] << "\t"
                    << bad_freq_ma[i] << "\t"
                    << bad_A1_freq[i] << "\t"
                    << bad_A2_freq[i] << "\t"
                    << bad_freq_value[i]
                    << endl;
        }
        obadsnp.close();
        LOGGER << bad_snp_freq.size() << " SNP(s) have large difference of allele frequency between the GWAS summary data and the reference sample. These SNPs have been saved in [" << badsnpfile << "]." << endl; 
    }

    if (_include.empty()) LOGGER.e(0, "none of the SNPs in the GWAS summary data can be found in the genotype data.");
    else LOGGER << _include.size() << " SNPs are matched to the genotype data." << endl;

    if (_mu.empty()) calcu_mu();

    _freq.resize(_include.size());
    _beta.resize(_include.size());
    _beta_se.resize(_include.size());
    _pval.resize(_include.size());
    _chisq.resize(_include.size());
    _N_o.resize(_include.size());
    for (i = 0; i < _include.size(); i++) {
        _freq[i] = freq_buf[indx[i]];
        _beta[i] = beta_buf[indx[i]];
        _beta_se[i] = beta_se_buf[indx[i]];
        chi_buf = _beta[i] / _beta_se[i];
        _chisq[i] = chi_buf * chi_buf;
        if (GC) {
            _pval[i] = StatFunc::pchisq(_chisq[i] / _GC_val, 1);
        } else _pval[i] = StatFunc::pchisq(_chisq[i], 1);
        _N_o[i] = N_o_buf[indx[i]];
    }
    _jma_Ve = _jma_Vp;
}

void gcta::init_massoc(string metafile, bool GC, double GC_val)
{
    read_metafile(metafile, GC, GC_val);

    int i = 0, j = 0, n = _keep.size(), m = _include.size();
    LOGGER << "Calculating the variance of SNP genotypes ..." << endl;
    _MSX_B.resize(m);
    _Nd.resize(m);

    if (_mu.empty()) calcu_mu();

    // Cache centered genotypes for COJO LD when the footprint is modest (~320MB cap).
    // Under N-way parallel fan-out, total RAM scales with the number of processes.
    const size_t cells = static_cast<size_t>(n) * static_cast<size_t>(m);
    const size_t max_cells = 40ull << 20; // ~40M doubles ≈ 320MB
    _cojo_X_ready = false;
    _cojo_X.resize(0, 0);
    _cojo_makex_cache_i = -1;
    _cojo_makex_cache.resize(0);
    _Z_cache_ready = false;
    _Z_cache.resize(0, 0);
    _Z_N_cache.resize(0, 0);
    _Z.resize(0, 0);
    _Z_N.resize(0, 0);
    if (cells > 0 && cells <= max_cells) {
        _cojo_X.resize(n, m);
        _cojo_X_ready = true;
    }

    #pragma omp parallel for
    for (i = 0; i < m; i++){
        eigenVector x;
        makex_eigenVector(i, x, true, true);
        _MSX_B[i] = x.squaredNorm() / (double)n;
        if (_cojo_X_ready) _cojo_X.col(i) = x;
    }
    if (_cojo_X_ready) {
        LOGGER << "COJO genotype cache enabled (" << m << " SNPs × " << n << " individuals)." << endl;
    }
    // Bit-packed keep-subset planes (primary LD kernel when float cache is off).
    build_cojo_bit_geno();
    if (_jma_actual_geno) {
        _MSX = _MSX_B;
        _Nd = _N_o;
    } else {
        _MSX = 2.0 * _freq.array()*(1.0 - _freq.array());
        for (i = 0; i < m; i++) _Nd[i] = (_jma_Vp - _MSX[i] * _beta[i] * _beta[i]) / (_MSX[i] * _beta_se[i] * _beta_se[i]) + 1; // revised by JY 25/11/13 according to Eq. 13, Yang et al. 2012 NG
    }
    build_cojo_bp_order();
}

void gcta::build_cojo_bit_geno()
{
    _cojo_bits_ready = false;
    _cojo_bit1.clear();
    _cojo_bit2.clear();
    _cojo_keep_mask.clear();
    _cojo_bit_mu.clear();
    _cojo_words = 0;
    const int n = static_cast<int>(_keep.size());
    const int m = static_cast<int>(_include.size());
    if (n < 1 || m < 1) return;

    _cojo_words = static_cast<size_t>((n + 63) / 64);
    _cojo_keep_mask.assign(_cojo_words, ~0ULL);
    if (n % 64 != 0) {
        _cojo_keep_mask.back() = (n % 64 == 0) ? ~0ULL : ((1ULL << (n % 64)) - 1ULL);
    }
    _cojo_bit1.resize(m);
    _cojo_bit2.resize(m);
    _cojo_bit_mu.resize(m);

    #pragma omp parallel for
    for (int j = 0; j < m; j++) {
        _cojo_bit1[j].assign(_cojo_words, 0ULL);
        _cojo_bit2[j].assign(_cojo_words, 0ULL);
        const int snp = _include[j];
        const double mu = _mu[snp];
        _cojo_bit_mu[j] = mu;
        const bool ref_is_a1 = (_allele1[snp] == _ref_A[snp]);
        for (int i = 0; i < n; i++) {
            const int indi = _keep[i];
            const bool s1 = _snp_1[snp][indi];
            const bool s2 = _snp_2[snp][indi];
            const bool missing = (s1 && !s2);
            const size_t w = static_cast<size_t>(i >> 6);
            const uint64_t bit = 1ULL << (i & 63);
            if (missing) {
                // Missing: leave bit1=0,bit2=0 here; handled via (~bit1|bit2) valid mask in dots.
                // Encode as bit1=1,bit2=0 to match PLINK missing so valid-mask excludes them.
                _cojo_bit1[j][w] |= bit;
            } else {
                // Dosage 0/1/2 in allele1 count, then flip if ref is allele2.
                int dose = (int)s1 + (int)s2;
                if (!ref_is_a1) dose = 2 - dose;
                // Store as: dose 2 → bit1=1,bit2=1; dose 1 → bit1=0,bit2=1; dose 0 → bit1=0,bit2=0
                if (dose >= 2) {
                    _cojo_bit1[j][w] |= bit;
                    _cojo_bit2[j][w] |= bit;
                } else if (dose == 1) {
                    _cojo_bit2[j][w] |= bit;
                }
            }
        }
    }
    _cojo_bits_ready = true;
    LOGGER << "COJO bit-packed genotypes enabled (" << m << " SNPs × " << n << " individuals, "
           << _cojo_words << " words/SNP)." << endl;
}

double gcta::cojo_snp_dot_bits(int i, int j) const
{
    const int n = static_cast<int>(_keep.size());
    const double mu_i = _cojo_bit_mu[i];
    const double mu_j = _cojo_bit_mu[j];
    int N_both = 0, S12 = 0, S1 = 0, S2 = 0;
    for (size_t w = 0; w < _cojo_words; w++) {
        const uint64_t mask = _cojo_keep_mask[w];
        if (mask == 0ULL) continue;
        uint64_t x11 = _cojo_bit1[i][w], x12 = _cojo_bit2[i][w];
        uint64_t x21 = _cojo_bit1[j][w], x22 = _cojo_bit2[j][w];
        // valid = not missing: missing is bit1=1,bit2=0 → (~bit1 | bit2) clears missing
        uint64_t both = (~x11 | x12) & (~x21 | x22) & mask;
        x11 &= both;
        x12 &= both;
        x21 &= both;
        x22 &= both;
        N_both += __builtin_popcountll(both);
        // dose = I(bit1)+I(bit2): dose2=(1,1), dose1=(0,1), dose0=(0,0)
        S12 += __builtin_popcountll(x11 & x21) + __builtin_popcountll(x11 & x22)
             + __builtin_popcountll(x12 & x21) + __builtin_popcountll(x12 & x22);
        S1 += __builtin_popcountll(x11) + __builtin_popcountll(x12);
        S2 += __builtin_popcountll(x21) + __builtin_popcountll(x22);
    }
    // GCTA makex: missing → centered 0; divide by full n (not N_both).
    return (S12 - mu_i * S2 - mu_j * S1 + (double)N_both * mu_i * mu_j) / (double)n;
}

double gcta::cojo_snp_dot(int i, int j)
{
    const int n = static_cast<int>(_keep.size());
    // Prefer bit kernel when available; float cache remains for MSX / small fixtures.
    if (_cojo_bits_ready) {
        return cojo_snp_dot_bits(i, j);
    }
    if (_cojo_X_ready) {
        return _cojo_X.col(i).dot(_cojo_X.col(j)) / (double)n;
    }
    // Cache-disabled fallback: hoist the outer SNP (first arg is fixed in init_B /
    // init_Z / insert_Z loops) so we do not rematerialize both vectors every pair.
    if (_cojo_makex_cache_i != i) {
        makex_eigenVector(i, _cojo_makex_cache, true, true);
        _cojo_makex_cache_i = i;
    }
    eigenVector x_j;
    makex_eigenVector(j, x_j, true, true);
    return _cojo_makex_cache.dot(x_j) / (double)n;
}

void gcta::build_cojo_bp_order()
{
    const int m = static_cast<int>(_include.size());
    _bp_order.resize(m);
    for (int i = 0; i < m; i++) _bp_order[i] = i;
    std::stable_sort(_bp_order.begin(), _bp_order.end(), [&](int a, int b) {
        const int ca = _chr[_include[a]], cb = _chr[_include[b]];
        if (ca != cb) return ca < cb;
        return _bp[_include[a]] < _bp[_include[b]];
    });
}

void gcta::cojo_window_bounds(int snp_include_idx, int &lo, int &hi) const
{
    // Half-open [lo, hi) into _bp_order for SNPs with abs(bp diff) < _jma_wind_size on same chr.
    const int m = static_cast<int>(_bp_order.size());
    lo = 0;
    hi = m;
    if (m == 0 || snp_include_idx < 0 || snp_include_idx >= (int)_include.size()) return;
    const int chr = _chr[_include[snp_include_idx]];
    const int bp = _bp[_include[snp_include_idx]];

    int c_lo = 0, c_hi = m;
    while (c_lo < m && _chr[_include[_bp_order[c_lo]]] < chr) c_lo++;
    c_hi = c_lo;
    while (c_hi < m && _chr[_include[_bp_order[c_hi]]] == chr) c_hi++;
    if (c_lo == c_hi) {
        lo = hi = c_lo;
        return;
    }

    // abs(bp_j - bp) < wind  <=>  bp - wind < bp_j < bp + wind
    const int left = bp - _jma_wind_size + 1; // first bp that can be included
    const int right = bp + _jma_wind_size - 1;
    auto bp_of = [&](int ord) { return _bp[_include[_bp_order[ord]]]; };

    int L = c_lo, R = c_hi;
    while (L < R) {
        int mid = (L + R) / 2;
        if (bp_of(mid) < left) L = mid + 1;
        else R = mid;
    }
    lo = L;
    L = c_lo;
    R = c_hi;
    while (L < R) {
        int mid = (L + R) / 2;
        if (bp_of(mid) <= right) L = mid + 1;
        else R = mid;
    }
    hi = L;
}

void gcta::read_fixed_snp(string snplistfile, string msg, vector<int> &pgiven, vector<int> &remain) {
    int i = 0, j = 0;
    vector<string> givenSNPs;
    read_snplist(snplistfile, givenSNPs, msg);
    if (givenSNPs.empty()) LOGGER.e(0, "failed to read any SNP from the file [" + snplistfile + "].");
    map<string, int> givenSNPs_map;
    pgiven.clear();
    remain.clear();
    for (i = 0; i < givenSNPs.size(); i++) givenSNPs_map.insert(pair<string, int>(givenSNPs[i], i));
    for (i = 0; i < _include.size(); i++) {
        if (givenSNPs_map.find(_snp_name[_include[i]]) != givenSNPs_map.end()) pgiven.push_back(i);
        else remain.push_back(i);
    }
    if (pgiven.size() > 0) LOGGER << pgiven.size() << " of them are matched to the genotype and summary data." << endl;
    else LOGGER.e(0, "none of the given SNPs can be matched to the genotype and summary data.");
}

void gcta::run_massoc_slct(string metafile, int wind_size, double p_cutoff, double collinear, int64_t top_SNPs, bool joint_only, bool GC, double GC_val, bool actual_geno, int mld_slct_alg)
{
    _jma_actual_geno = actual_geno;
    _jma_wind_size = wind_size;
    _jma_p_cutoff = p_cutoff;
    _jma_collinear = collinear;
    _jma_snpnum_backward = 0;
    _jma_snpnum_collienar = 0;
    init_massoc(metafile, GC, GC_val);
    if (top_SNPs < 0) top_SNPs = 1e10;
    else {
        _jma_p_cutoff = 0.5;
        LOGGER << "The threshold p-value has been set to 0.5 because of the --cojo-top-SNPs option." << endl;
    }
    int i = 0, j = 0;
    vector<int> slct, remain;
    eigenVector bC, bC_se, pC;
    LOGGER << endl;
    if (!joint_only && mld_slct_alg<2) {
        LOGGER << "Performing "<< ((mld_slct_alg==0)?"stepwise":"forward") << " model selection on " << _include.size() << " SNPs to select association signals ... (p cutoff = " << _jma_p_cutoff << "; ";
        LOGGER << "collinearity cutoff = " << _jma_collinear << ")"<< endl;
        if (!_jma_actual_geno) LOGGER << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;
        stepwise_slct(slct, remain, bC, bC_se, pC, mld_slct_alg, top_SNPs);
        if (slct.empty()) {
            LOGGER << "No SNPs have been selected." << endl;
            return;
        }
    } else { // mld_slct_alg = 2 for backward selection
        for (i = 0; i < _include.size(); i++) slct.push_back(i);
        if (mld_slct_alg==2) {
            LOGGER << "Performing backward selection on " << _include.size() << " SNPs at threshold p-value = " << _jma_p_cutoff << " ..." << endl;
            slct_stay(slct, bC, bC_se, pC);
        }
    }

    // joint analysis
    eigenVector bJ, bJ_se, pJ;
    LOGGER << "Performing joint analysis on all the " << slct.size();
    if (joint_only) LOGGER << " SNPs ..." << endl;
    else LOGGER << " selected signals ..." << endl;
    if (slct.size() >= _keep.size()) LOGGER.e(0, "too many SNPs. The number of SNPs in a joint analysis should not be larger than the sample size.");
    massoc_joint(slct, bJ, bJ_se, pJ);
    eigenMatrix rval(slct.size(), slct.size());
    LD_rval(slct, rval);
    if (_jma_actual_geno) LOGGER << "Residual variance = " << _jma_Ve << endl;
    massoc_slct_output(joint_only, slct, bJ, bJ_se, pJ, rval);

    // output conditional results
    if (!joint_only && mld_slct_alg!=2) {
        massoc_cond_output(remain, bC, bC_se, pC);
        LOGGER << "(" << _jma_snpnum_backward << " SNPs eliminated by backward selection and " << _jma_snpnum_collienar << " SNPs filtered by collinearity test are not included in the output)" << endl;
    }
}

void gcta::run_massoc_cond(string metafile, string snplistfile, int wind_size, double collinear, bool GC, double GC_val, bool acutal_geno) {
    _jma_actual_geno = acutal_geno;
    _jma_wind_size = wind_size;
    _jma_collinear = collinear;
    init_massoc(metafile, GC, GC_val);
    vector<int> pgiven, remain;
    read_fixed_snp(snplistfile, "given SNPs", pgiven, remain);

    eigenVector bC, bC_se, pC;
    LOGGER << "Performing single-SNP association analysis conditional on the " << pgiven.size() << " given SNPs ... ";
    LOGGER << "(collinearity cutoff = " << _jma_collinear << ")" << endl;
    if (!_jma_actual_geno) LOGGER << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;

    string filename = _out + ".given.cojo";
    LOGGER<<"Saving the summary statistics of the given SNPs in the file [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << endl;
    int i = 0, j = 0;
    for (i = 0; i < pgiven.size(); i++) {
        j = pgiven[i];
        ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t" << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t" << _pval[j] << endl;
    }
    ofile.close();
    
    massoc_cond(pgiven, remain, bC, bC_se, pC);
    massoc_cond_output(remain, bC, bC_se, pC);
}

void gcta::massoc_slct_output(bool joint_only, vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ, eigenMatrix &rval)
{
    string filename = _out + ".jma.cojo";
    if (joint_only) LOGGER << "Saving the joint analysis result of " << slct.size() << " SNPs to [" + filename + "] ..." << endl;
    else LOGGER << "Saving the " << slct.size() << " independent signals to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << ((_GC_val > 0) ? "_GC" : "") << "\tn\tfreq_geno\tbJ\tbJ_se\tpJ" << ((_GC_val > 0) ? "_GC" : "") << "\tLD_r" << endl;
    int i = 0, j = 0;
    for (i = 0; i < slct.size(); i++) {
        j = slct[i];
        ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t";
        ofile << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t";
        ofile << _pval[j] << "\t" << _Nd[j] << "\t" << 0.5 * _mu[_include[j]] << "\t" << bJ[i] << "\t" << bJ_se[i] << "\t" << pJ[i] << "\t";
        if (i == slct.size() - 1) ofile << 0 << endl;
        else ofile << rval(i, i + 1) << endl;
    }
    ofile.close();

    filename = _out + ".ldr.cojo";
    if (joint_only) LOGGER << "Saving the LD structure of " << slct.size() << " SNPs to [" + filename + "] ..." << endl;
    else LOGGER << "Saving the LD structure of " << slct.size() << " independent signals to [" + filename + "] ..." << endl;
    ofile.open(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "SNP\t";
    for (i = 0; i < slct.size(); i++) ofile << _snp_name[_include[slct[i]]] << "\t";
    ofile << endl;
    for (i = 0; i < slct.size(); i++) {
        ofile << _snp_name[_include[slct[i]]] << "\t";
        for (j = 0; j < slct.size(); j++) ofile << rval(i, j) << "\t";
        ofile << endl;
    }
    ofile.close();
}

void gcta::set_massoc_pC_thresh(double thresh){
    g_massoc_out_thresh = thresh;
}

void gcta::massoc_cond_output(vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC)
{
    double out_thresh = g_massoc_out_thresh;
    int i = 0, j = 0;
    string filename = _out + ".cma.cojo";
    LOGGER << "Saving the conditional analysis results of " << remain.size() << " remaining SNPs to [" + filename + "] ..." << endl;
    ofstream ofile(filename.c_str());
    if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
    ofile << "Chr\tSNP\tbp\trefA\tfreq\tb\tse\tp" << ((_GC_val > 0) ? "_GC" : "") << "\tn\tfreq_geno\tbC\tbC_se\tpC" << ((_GC_val > 0) ? "_GC" : "") << endl;

    if(out_thresh < 0){
        for (i = 0; i < remain.size(); i++) {
            j = remain[i];
            ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t";
            ofile << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t";
            ofile << _pval[j] << "\t" << _Nd[j] << "\t" << 0.5 * _mu[_include[j]] << "\t";
            if (pC[i] > 1.5) ofile << "NA\tNA\tNA" << endl;
            else ofile << bC[i] << "\t" << bC_se[i] << "\t" << pC[i] << endl;
        }
    }else{
        LOGGER << "Restricting output threshold to " << out_thresh << "." << endl;
        for (i = 0; i < remain.size(); i++) {
            if(pC[i] < out_thresh){
                j = remain[i];
                ofile << _chr[_include[j]] << "\t" << _snp_name[_include[j]] << "\t" << _bp[_include[j]] << "\t";
                ofile << _ref_A[_include[j]] << "\t" << _freq[j] << "\t" << _beta[j] << "\t" << _beta_se[j] << "\t";
                ofile << _pval[j] << "\t" << _Nd[j] << "\t" << 0.5 * _mu[_include[j]] << "\t";
                ofile << bC[i] << "\t" << bC_se[i] << "\t" << pC[i] << endl;
            }
        }
    }
    ofile.close();
}

void gcta::stepwise_slct(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC, int mld_slct_alg, uint64_t top_SNPs)
{
    int i = 0, i_buf = 0;
    vector<double> p_buf, chisq;
    
    for (i = 0; i < _include.size(); i++) {
        remain.push_back(i);
    }
    
    for (i = 0; i < remain.size(); i++) {
        chisq.push_back(_chisq[remain[i]]);
    }

    int max_chisq_index = max_element(chisq.begin(), chisq.end()) - chisq.begin();
    if (_pval[remain[max_chisq_index]] >= _jma_p_cutoff) {
        return;
    } else {
        slct.push_back(remain[max_chisq_index]);
    }

    remain.erase(remain.begin() + max_chisq_index);
    
    int prev_num = 0;
    if (mld_slct_alg==0 && _jma_p_cutoff > 1e-3){
        LOGGER << "Switched to perform a forward model selection because the significance level is too low..." << endl;
        mld_slct_alg=1;
    }
    bool slct_only_contain_remain = true;
    while (!remain.empty()) {
        if (slct.size() > 0) {
            slct_only_contain_remain = false;
            if (slct_entry(slct, remain, bC, bC_se, pC)) {
                if (mld_slct_alg == 0) slct_stay(slct, bC, bC_se, pC);
            } else {
                break;
            }
        } else {
            chisq.clear();
            for (i = 0; i < remain.size(); i++) {
                chisq.push_back(_chisq[remain[i]]);
            }
            int max_chisq_index =
                max_element(chisq.begin(), chisq.end()) - chisq.begin();
            if (_pval[remain[max_chisq_index]] >= _jma_p_cutoff) {
                break;
            } else {
                slct.push_back(remain[max_chisq_index]);
            }
            remain.erase(remain.begin() + max_chisq_index);
            slct_only_contain_remain = true;
        }
        if (!slct_only_contain_remain) {
            if (slct.size() % 5 == 0 && slct.size() > prev_num)
                LOGGER << slct.size() << " associated SNPs have been selected."
                       << endl;
            if (slct.size() > prev_num) prev_num = slct.size();
            if (slct.size() >= top_SNPs) break;
        }
    }
    if (slct_only_contain_remain) {
        if (mld_slct_alg == 0) {
            slct_stay(slct, bC, bC_se, pC);
        } else {
            LOGGER.e(0, "Warning, this should not happend");
        }
    }

    if (_jma_p_cutoff > 1e-3) {
        LOGGER << "Performing backward elimination..." << endl;
        slct_stay(slct, bC, bC_se, pC);
    }
    LOGGER << "Finally, " << slct.size() << " associated SNPs are selected." << endl;
}

bool gcta::slct_entry(vector<int> &slct, vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC) {
    if (slct.size() < 1) {
        LOGGER.e(0, "Warning, this should not happend");
    } 
    int i = 0, m = 0;
    massoc_cond(slct, remain, bC, bC_se, pC);
    vector<double> pC_buf;
    eigenVector2Vector(pC, pC_buf);
    while (true) {
        m = min_element(pC_buf.begin(), pC_buf.end()) - pC_buf.begin();
        if (pC_buf[m] >= _jma_p_cutoff) return (false);
        if (insert_B_and_Z(slct, remain[m])) {
            slct.push_back(remain[m]);
            stable_sort(slct.begin(), slct.end());
            remain.erase(remain.begin() + m);
            return (true);
        }
        pC_buf.erase(pC_buf.begin() + m);
        remain.erase(remain.begin() + m);
    }
}

void gcta::slct_stay(vector<int> &slct, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ) {
    if (_B_N.cols() < 1) {
        if (!init_B(slct)) LOGGER.e(0, "there is a collinearity problem of the given list of SNPs.\nYou can try the option --cojo-slct to remove one of each pair of highly correlated SNPs.");
    }

    vector<double> pJ_buf;
    while (!slct.empty()) {
        massoc_joint(slct, bJ, bJ_se, pJ);
        eigenVector2Vector(pJ, pJ_buf);
        int m = max_element(pJ_buf.begin(), pJ_buf.end()) - pJ_buf.begin();
        if (pJ[m] > _jma_p_cutoff) {
            _jma_snpnum_backward++;
            erase_B_and_Z(slct, slct[m]);
            slct.erase(slct.begin() + m);
        } else break;
    }
}

void gcta::eigenVector2Vector(eigenVector &x, vector<double> &y) {
    y.resize(x.size());
    for (int i = 0; i < x.size(); i++) y[i] = x[i];
}

double gcta::massoc_calcu_Ve(const vector<int> &slct, eigenVector &bJ, eigenVector &b) {
    double Ve = 0.0;
    int n = bJ.size();
    vector<double> Nd_buf(n);
    for (int k = 0; k < n; k++) {
        Nd_buf[k] = _Nd[slct[k]];
        Ve += _D_N[k] * bJ[k] * b[k];
    }
    double d_buf = CommFunc::median(Nd_buf);
    if (d_buf - n < 1) LOGGER.e(0, "no degree of freedom is left for the residues. The model is over-fitted. Please specify a more stringent p-value cut-off.");
    Ve = ((d_buf - 1) * _jma_Vp - Ve) / (d_buf - n);
    if (Ve <= 0.0) LOGGER.e(0, "residual variance is out of boundary. The model is over-fitted. Please specify a more stringent p-value cut-off.");
    return Ve;
}

void gcta::massoc_joint(const vector<int> &indx, eigenVector &bJ, eigenVector &bJ_se, eigenVector &pJ) {
    if (_B_N.cols() < 1) {
        if (!init_B(indx)) LOGGER.e(0, "there is a collinearity problem of the given list of SNPs.\nYou can try the option --cojo-slct to remove one of each pair of highly correlated SNPs.");
    }

    int i = 0, n = indx.size();
    double chisq = 0.0;
    eigenVector b(n);
    for (i = 0; i < n; i++) b[i] = _beta[indx[i]];
    bJ.resize(n);
    bJ_se.resize(n);
    pJ.resize(n);
    bJ = _B_N_i * _D_N.asDiagonal() * b;
    bJ_se = _B_N_i.diagonal();
    pJ = eigenVector::Ones(n);
    if (_jma_actual_geno) _jma_Ve = massoc_calcu_Ve(indx, bJ, b);
    bJ_se *= _jma_Ve;
    for (i = 0; i < n; i++) {
        if (bJ_se[i] > 1.0e-30) {
            bJ_se[i] = sqrt(bJ_se[i]);
            chisq = bJ[i] / bJ_se[i];
            if (_GC_val > 0) pJ[i] = StatFunc::pchisq(chisq * chisq / _GC_val, 1);
            else pJ[i] = StatFunc::pchisq(chisq*chisq, 1);
        } else {
            bJ[i] = 0.0;
            bJ_se[i] = 0.0;
        }
    }
}

void gcta::massoc_cond(const vector<int> &slct, const vector<int> &remain, eigenVector &bC, eigenVector &bC_se, eigenVector &pC) {
    if (slct.size() < 1) {
        LOGGER.e(0, "Warning, this should not happend");
    }
    if (_B_N.cols() < 1) {
        if (!init_B(slct)) LOGGER.e(0, "there is a collinearity problem of the given list of SNPs.\nYou can try the option --cojo-slct to remove one of each pair of highly correlated SNPs.");
    }
    if (!_Z_cache_ready) init_Z(slct);

    int i = 0, j = 0, n = slct.size();
    double chisq = 0.0;
    eigenVector b(n), se(n);
    for (i = 0; i < n; i++) {
        b[i] = _beta[slct[i]];
        se[i] = _beta_se[slct[i]];
    }
    eigenVector bJ1 = _B_N_i * _D_N.asDiagonal() * b;
    if (_jma_actual_geno) _jma_Ve = massoc_calcu_Ve(slct, bJ1, b);

    double B2 = 0.0;
    bC = eigenVector::Zero(remain.size());
    bC_se = eigenVector::Zero(remain.size());
    pC = eigenVector::Constant(remain.size(), 2);
    eigenVector Z_Bi(n), Z_Bi_buf(n);
    double cutoff = 1e-10 * _jma_Vp;
    //LOGGER << "Cutoff of bC_se " << cutoff << endl;
    for (i = 0; i < remain.size(); i++) {
        j = remain[i];
        B2 = _MSX[j] * _Nd[j];
        if (!CommFunc::FloatEqual(B2, 0.0)) {
            Z_Bi = _Z_N_cache.col(j).transpose() * _B_N_i;
            Z_Bi_buf = _Z_cache.col(j).transpose() * _B_i;
            if (_Z_cache.col(j).dot(Z_Bi_buf) / _MSX_B[j] < _jma_collinear) {
                bC[i] = _beta[j] - Z_Bi.cwiseProduct(_D_N).dot(b) / B2;
                bC_se[i] = 1/B2;     // Revised by Zhihong 4 April 2017 //bC_se[i] = (B2 - _Z_N.col(j).dot(Z_Bi)) / (B2 * B2);
            }
        }
        if (_jma_actual_geno) bC_se[i] *= _jma_Ve - (B2 * bC[i] * _beta[j]) / (_Nd[j] - n - 1);
        else bC_se[i] *= _jma_Ve;
        if (bC_se[i] > cutoff) {
            bC_se[i] = sqrt(bC_se[i]);
            chisq = bC[i] / bC_se[i];
            if (_GC_val > 0) pC[i] = StatFunc::pchisq(chisq * chisq / _GC_val, 1);
            else pC[i] = StatFunc::pchisq(chisq*chisq, 1);
        }
    }
}

bool gcta::init_B(const vector<int> &indx)
{
    if (indx.size() < 1) {
        LOGGER.e(0, "slct size is zero will cause Eigen Matrix of"
            "Vector operation error.");
    }
    int i = 0, j = 0, k = 0, n = _keep.size();
    double d_buf = 0.0;
    eigenVector diagB(indx.size());
    _B.resize(indx.size(), indx.size());
    _B_N.resize(indx.size(), indx.size());
    _D_N.resize(indx.size());
    for (i = 0; i < indx.size(); i++) {
        _D_N[i] = _MSX[indx[i]] * _Nd[indx[i]];
        _B.startVec(i);
        _B_N.startVec(i);
        _B.insertBack(i, i) = _MSX_B[indx[i]];
        _B_N.insertBack(i, i) = _D_N[i];
        diagB[i] = _MSX_B[indx[i]];
        for (j = i + 1; j < indx.size(); j++) {
            if (_jma_actual_geno || (_chr[_include[indx[i]]] == _chr[_include[indx[j]]] && abs(_bp[_include[indx[i]]] - _bp[_include[indx[j]]]) < _jma_wind_size)) {
                d_buf = cojo_snp_dot(indx[i], indx[j]);
                _B.insertBack(j, i) = d_buf;
                _B_N.insertBack(j, i) = d_buf * min(_Nd[indx[i]], _Nd[indx[j]]) * sqrt(_MSX[indx[i]] * _MSX[indx[j]] / (_MSX_B[indx[i]] * _MSX_B[indx[j]]));
            }
        }
    }
    _B.finalize();
    _B_N.finalize();

    SimplicialLDLT<eigenSparseMat> ldlt_B(_B);

    if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30) return false;

    _B_i.resize(indx.size(), indx.size());
    _B_i.setIdentity();
    _B_i = ldlt_B.solve(_B_i).eval();
    if ((1 - eigenVector::Constant(indx.size(), 1).array() / (diagB.array() * _B_i.diagonal().array())).maxCoeff() > _jma_collinear) return false;
    SimplicialLDLT<eigenSparseMat> ldlt_B_N(_B_N);
    _B_N_i.resize(indx.size(), indx.size());
    _B_N_i.setIdentity();
    _B_N_i = ldlt_B_N.solve(_B_N_i).eval();
    return true;
}

void gcta::fill_z_row_for_selected(int row, int selected_include_idx)
{
    const int m = static_cast<int>(_include.size());
    if (row < 0 || row >= _Z_cache.rows() || selected_include_idx < 0 || selected_include_idx >= m) return;

    int lo = 0, hi = m;
    if (_jma_actual_geno) {
        lo = 0;
        hi = static_cast<int>(_bp_order.size());
    } else {
        cojo_window_bounds(selected_include_idx, lo, hi);
    }

    // Thread-safe: bits / float cache are read-only; avoid shared makex hoist under OpenMP.
    #pragma omp parallel for schedule(static)
    for (int t = lo; t < hi; t++) {
        const int j = _bp_order[t];
        if (!_jma_actual_geno) {
            if (selected_include_idx == j) continue;
            if (_chr[_include[selected_include_idx]] != _chr[_include[j]]) continue;
            if (abs(_bp[_include[selected_include_idx]] - _bp[_include[j]]) >= _jma_wind_size) continue;
        }
        const double d_buf = cojo_snp_dot(j, selected_include_idx);
        _Z_cache(row, j) = d_buf;
        _Z_N_cache(row, j) = d_buf * min(_Nd[selected_include_idx], _Nd[j])
            * sqrt(_MSX[selected_include_idx] * _MSX[j] / (_MSX_B[selected_include_idx] * _MSX_B[j]));
    }
}

void gcta::append_z_cache_row(int insert_include_idx, int sorted_pos)
{
    const int k = static_cast<int>(_Z_cache.rows());
    const int m = static_cast<int>(_include.size());
    _Z_cache.conservativeResize(k + 1, m);
    _Z_N_cache.conservativeResize(k + 1, m);
    _Z_cache.row(k).setZero();
    _Z_N_cache.row(k).setZero();
    fill_z_row_for_selected(k, insert_include_idx);
    if (sorted_pos != k) {
        // Move last row to sorted_pos (same append-then-permute pattern as B).
        const eigenVector z_last = _Z_cache.row(k);
        const eigenVector zn_last = _Z_N_cache.row(k);
        for (int r = k; r > sorted_pos; r--) {
            _Z_cache.row(r) = _Z_cache.row(r - 1);
            _Z_N_cache.row(r) = _Z_N_cache.row(r - 1);
        }
        _Z_cache.row(sorted_pos) = z_last;
        _Z_N_cache.row(sorted_pos) = zn_last;
    }
}

void gcta::drop_z_cache_row(int erase_pos)
{
    const int k = static_cast<int>(_Z_cache.rows());
    const int m = static_cast<int>(_Z_cache.cols());
    if (k <= 1) {
        _Z_cache.resize(0, 0);
        _Z_N_cache.resize(0, 0);
        _Z_cache_ready = false;
        return;
    }
    eigenMatrix Z_new(k - 1, m), ZN_new(k - 1, m);
    if (erase_pos > 0) {
        Z_new.topRows(erase_pos) = _Z_cache.topRows(erase_pos);
        ZN_new.topRows(erase_pos) = _Z_N_cache.topRows(erase_pos);
    }
    if (erase_pos < k - 1) {
        Z_new.bottomRows(k - 1 - erase_pos) = _Z_cache.bottomRows(k - 1 - erase_pos);
        ZN_new.bottomRows(k - 1 - erase_pos) = _Z_N_cache.bottomRows(k - 1 - erase_pos);
    }
    _Z_cache.swap(Z_new);
    _Z_N_cache.swap(ZN_new);
}

void gcta::init_Z(const vector<int> &indx)
{
    if (indx.size() < 1) {
        LOGGER.e(0,
                 "slct size is zero will cause Eigen Matrix of Vector "
                 "operation error.");
    }
    const int k = static_cast<int>(indx.size());
    const int m = static_cast<int>(_include.size());
    _Z_cache = eigenMatrix::Zero(k, m);
    _Z_N_cache = eigenMatrix::Zero(k, m);
    for (int i = 0; i < k; i++) {
        fill_z_row_for_selected(i, indx[i]);
    }
    _Z_cache_ready = true;
    // Keep sparse _Z/_Z_N empty; massoc_cond reads dense caches.
    _Z.resize(0, 0);
    _Z_N.resize(0, 0);
}

bool gcta::insert_B_and_Z(const vector<int> &indx, int insert_indx)
{
    if (indx.size() < 1) {
        LOGGER.e(0,
                 "slct size is zero will cause Eigen Matrix of Vector "
                 "operation error.");
    }
    int i = 0, j = 0, n = _keep.size();
    double d_buf = 0.0;
    vector<int> ix(indx);
    ix.push_back(insert_indx);
    stable_sort(ix.begin(), ix.end());
    eigenSparseMat B_buf(_B), B_N_buf(_B_N);
    _B.resize(ix.size(), ix.size());
    _B_N.resize(ix.size(), ix.size());
    bool get_insert_col = false, get_insert_row = false;
    int pos = find(ix.begin(), ix.end(), insert_indx) - ix.begin();
    eigenVector diagB(ix.size());
    for (j = 0; j < ix.size(); j++) {
        _B.startVec(j);
        _B_N.startVec(j);
        _B.insertBack(j, j) = _MSX_B[ix[j]];
        _B_N.insertBack(j, j) = _MSX[ix[j]] * _Nd[ix[j]];
        diagB[j] = _MSX_B[ix[j]];
        if (insert_indx == ix[j]) get_insert_col = true;
        get_insert_row = get_insert_col;
        for (i = j + 1; i < ix.size(); i++) {
            if (insert_indx == ix[i]) get_insert_row = true;
            if (insert_indx == ix[j] || insert_indx == ix[i]) {
                if (_jma_actual_geno || (_chr[_include[ix[i]]] == _chr[_include[ix[j]]] && abs(_bp[_include[ix[i]]] - _bp[_include[ix[j]]]) < _jma_wind_size)) {
                    d_buf = cojo_snp_dot(ix[i], ix[j]);
                    _B.insertBack(i, j) = d_buf;
                    _B_N.insertBack(i, j) = d_buf * min(_Nd[ix[i]], _Nd[ix[j]]) * sqrt(_MSX[ix[i]] * _MSX[ix[j]] / (_MSX_B[ix[i]] * _MSX_B[ix[j]]));
                }
            } else {
                if (B_buf.coeff(i - get_insert_row, j - get_insert_col) != 0) {
                    _B.insertBack(i, j) = B_buf.coeff(i - get_insert_row, j - get_insert_col);
                    _B_N.insertBack(i, j) = B_N_buf.coeff(i - get_insert_row, j - get_insert_col);
                }
            }
        }
    }
    _B.finalize();
    _B_N.finalize();

    // Rank-1 inverse update when a prior inverse matches the pre-insert selection.
    // Fall back to LDLT when sizes don't match (first build after init_B edge cases).
    const bool can_incr = (_B_i.rows() == (int)indx.size() && _B_i.cols() == (int)indx.size()
                           && _B_N_i.rows() == (int)indx.size() && _B_N_i.cols() == (int)indx.size()
                           && indx.size() > 0);
    eigenMatrix B_i_buf = _B_i;
    eigenMatrix B_N_i_buf = _B_N_i;

    if (can_incr) {
        const int n_new = static_cast<int>(ix.size());
        eigenVector cB = extract_sparse_col_excl_diag(_B, pos, n_new);
        eigenVector cBN = extract_sparse_col_excl_diag(_B_N, pos, n_new);
        const double dB = _MSX_B[insert_indx];
        const double dBN = _MSX[insert_indx] * _Nd[insert_indx];

        eigenMatrix B_i_app, B_N_i_app;
        // Reject non-PD / singular Schur complement explicitly (stock LDLT minCoeff < 0).
        if (!inv_update_forward_append(B_i_buf, cB, dB, B_i_app)
            || !inv_update_forward_append(B_N_i_buf, cBN, dBN, B_N_i_app)) {
            _jma_snpnum_collienar++;
            _B = B_buf;
            _B_N = B_N_buf;
            _B_i = B_i_buf;
            _B_N_i = B_N_i_buf;
            return false;
        }
        _B_i = move_last_to_pos(B_i_app, pos);
        _B_N_i = move_last_to_pos(B_N_i_app, pos);

        // Same three guards as the LDLT path: PD, cond(D) > 30, per-SNP collinearity.
        SimplicialLDLT<eigenSparseMat> ldlt_B(_B);
        if (ldlt_B.vectorD().minCoeff() < 0
            || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30
            || (1 - eigenVector::Constant(ix.size(), 1).array() / (diagB.array() * _B_i.diagonal().array())).maxCoeff() > _jma_collinear) {
            _jma_snpnum_collienar++;
            _B = B_buf;
            _B_N = B_N_buf;
            _B_i = B_i_buf;
            _B_N_i = B_N_i_buf;
            return false;
        }
    } else {
        SimplicialLDLT<eigenSparseMat> ldlt_B(_B);
        _B_i.resize(ix.size(), ix.size());
        _B_i.setIdentity();
        _B_i = ldlt_B.solve(_B_i).eval();
        if (ldlt_B.vectorD().minCoeff() < 0 || sqrt(ldlt_B.vectorD().maxCoeff() / ldlt_B.vectorD().minCoeff()) > 30 || (1 - eigenVector::Constant(ix.size(), 1).array() / (diagB.array() * _B_i.diagonal().array())).maxCoeff() > _jma_collinear) {
            _jma_snpnum_collienar++;
            _B = B_buf;
            _B_N = B_N_buf;
            _B_i = B_i_buf;
            _B_N_i = B_N_i_buf;
            return false;
        }
        SimplicialLDLT<eigenSparseMat> ldlt_B_N(_B_N);
        _B_N_i.resize(ix.size(), ix.size());
        _B_N_i.setIdentity();
        _B_N_i = ldlt_B_N.solve(_B_N_i).eval();
    }
    _D_N.resize(ix.size());
    for (j = 0; j < ix.size(); j++) {
        _D_N[j] = _MSX[ix[j]] * _Nd[ix[j]];
    }

    // Incremental dense Z row append (OpenMP window fill); skip if Z never built.
    if (!_Z_cache_ready) return true;
    append_z_cache_row(insert_indx, pos);

    return true;
}

void gcta::erase_B_and_Z(const vector<int> &indx, int erase_indx) {
    int i = 0, j = 0;   
    eigenSparseMat B_dense(_B), B_N_dense(_B_N);
    _B.resize(indx.size() - 1, indx.size() - 1);
    _B_N.resize(indx.size() - 1, indx.size() - 1);
    _D_N.resize(indx.size() - 1);
    int pos = find(indx.begin(), indx.end(), erase_indx) - indx.begin();
    bool get_insert_col = false, get_insert_row = false;
    for (j = 0; j < indx.size(); j++) {
        if (erase_indx == indx[j]) {
            get_insert_col = true;
            continue;
        }
        _B.startVec(j - get_insert_col);
        _B_N.startVec(j - get_insert_col);
        _D_N[j - get_insert_col] = _MSX[indx[j]] * _Nd[indx[j]];
        get_insert_row = get_insert_col;
        for (i = j; i < indx.size(); i++) {
            if (erase_indx == indx[i]) {
                get_insert_row = true;
                continue;
            }
            if (B_dense.coeff(i, j) != 0) {
                _B.insertBack(i - get_insert_row, j - get_insert_col) = B_dense.coeff(i, j);
                _B_N.insertBack(i - get_insert_row, j - get_insert_col) = B_N_dense.coeff(i, j);
            }
        }
    }
    _B.finalize();
    _B_N.finalize();

    // Rank-1 downdate of inverses when sizes match; else LDLT rebuild.
    // Always refresh inverses (even if _Z was never built) so backward elimination stays consistent.
    const bool can_decr = (_B_i.rows() == (int)indx.size() && _B_i.cols() == (int)indx.size()
                           && _B_N_i.rows() == (int)indx.size() && _B_N_i.cols() == (int)indx.size()
                           && indx.size() > 0);
    if (can_decr) {
        eigenMatrix B_i_new, B_N_i_new;
        inv_update_backward(_B_i, pos, B_i_new);
        inv_update_backward(_B_N_i, pos, B_N_i_new);
        _B_i.swap(B_i_new);
        _B_N_i.swap(B_N_i_new);
    } else {
        _B_i.resize(indx.size() - 1, indx.size() - 1);
        _B_i.setIdentity();
        _B_N_i.resize(indx.size() - 1, indx.size() - 1);
        _B_N_i.setIdentity();
        if (indx.size() > 1) {
            SimplicialLDLT<eigenSparseMat> ldlt_B(_B);
            _B_i = ldlt_B.solve(_B_i).eval();
            SimplicialLDLT<eigenSparseMat> ldlt_B_N(_B_N);
            _B_N_i = ldlt_B_N.solve(_B_N_i).eval();
        }
    }

    if (!_Z_cache_ready) return;
    drop_z_cache_row(pos);
}

/*
double gcta::crossprod(vector<float> &x_i, vector<float> &x_j) {
    double prod = 0.0;
    for (int i = 0; i < x_i.size(); i++) prod += x_i[i] * x_j[i];
    return (prod / x_i.size());
}
*/

void gcta::LD_rval(const vector<int> &indx, eigenMatrix &rval) {
    int i = 0, j = 0;
    eigenVector sd(indx.size());
    for (i = 0; i < indx.size(); i++) sd[i] = sqrt(_MSX_B[indx[i]]);
    for (j = 0; j < indx.size(); j++) {
        rval(j, j) = 1.0;
        for (i = j + 1; i < indx.size(); i++) rval(i, j) = rval(j, i) = _B.coeff(i, j) / sd[i] / sd[j];
    }
}

void gcta::run_massoc_sblup(string metafile, int wind_size, double lambda)
{
    _jma_wind_size = wind_size;
    init_massoc(metafile, false, -1);

    int j = 0;
    eigenVector bJ;
    LOGGER << "\nPerforming joint analysis on all the " << _include.size() << " SNPs ..." << endl;
    LOGGER << "(Assuming complete linkage equilibrium between SNPs which are more than " << _jma_wind_size / 1e6 << "Mb away from each other)" << endl;
    if (massoc_sblup(lambda, bJ)) {
        string filename = _out + ".sblup.cojo";
        LOGGER << "Saving the joint analysis result of " << _include.size() << " SNPs to [" + filename + "] ..." << endl;
        ofstream ofile(filename.c_str());
        if (!ofile) LOGGER.e(0, "cannot open the file [" + filename + "] to write.");
        for (j = 0; j < _include.size(); j++) {
            // change here: convert from u to b before writing to file
            ofile << _snp_name[_include[j]] << "\t" << _ref_A[_include[j]] << "\t" << _beta[j] << "\t" << bJ[j] / sqrt(_MSX[j]) << endl;
        }
        ofile.close();
    } else LOGGER.e(0, "Jacobi iteration cannot converge. You can increase the maximum number of iterations by the option --massoc-sblup-maxit.");
}

bool gcta::massoc_sblup(double lambda, eigenVector &bJ)
{
    int i = 0, j = 0, k = 0, n = _keep.size(), m = _include.size();
    //eigenVector D(_include.size());
    //eigenVector invsqrtvar(_include.size());
    eigenVector Xty(_include.size());
    uint64_t num_err_snp = 0;
    for(int i = 0; i < _include.size(); i++){
        //D[i] = _MSX[i] * _N_o[i];
        //invsqrtvar[i] = 1/sqrt(_MSX[i]);
        Xty[i] = sqrt(_MSX[i]) * _N_o[i] * _beta[i];
        if(_MSX_B[i] < 1e-10){
            num_err_snp++;
        }
    }
    if(num_err_snp > 0){
        LOGGER.e(0, "there are " + to_string(num_err_snp) + " SNPs with MAF=0.");
    }

    //eigenVector Xty = invsqrtvar.array() * D.array() * _beta.array();


    /*
    vector<int> nz_i(m), nz_j(m);
    for (i = 0; i < _include.size(); i++){
        nz_i[i]=0;
        nz_j[i]=0;
    }
    */
    LOGGER << "Calculating the LD correlation matrix of all the " << _include.size() << " SNPs..." << endl;
    eigenSparseMat B(_include.size(), _include.size());
    vector<int> n_elements(_include.size());

    for (i = 0; i < _include.size(); i++) {
        int count = 1; 
        for (j = i + 1; j < _include.size(); j++) {
            if (_chr[_include[i]] == _chr[_include[j]] && abs(_bp[_include[i]] - _bp[_include[j]]) < _jma_wind_size) {
                count++;
            }
        }
        n_elements[i] = count;
    }
    B.reserve(n_elements);

    //ofstream out_debug(_out + ".debug");
    for (i = 0; i < _include.size(); i++) {
        // change here: get standardized genotypes
        eigenVector x_i(_keep.size());
        makex_eigenVector_std(i, x_i, false, sqrt(_MSX_B[i]));
        B.insertBackUncompressed(i, i) =  _N_o[i] + lambda;
       // out_debug << i << "\t" << i << "\t" << _N_o[i] + lambda << endl;
        #pragma omp parallel for ordered schedule(static) private(j)
        for (j = i + 1; j < _include.size(); j++) {
            if (_chr[_include[i]] == _chr[_include[j]] && abs(_bp[_include[i]] - _bp[_include[j]]) < _jma_wind_size) {
                eigenVector x_j(_keep.size());
                // change here: get standardized genotypes
                makex_eigenVector_std(j, x_j, false, sqrt(_MSX_B[j]));
                double prod = x_i.dot(x_j) / (double)n;
                double temp = prod * min(_Nd[i], _Nd[j]) * sqrt(_MSX[i] * _MSX[j] / (_MSX_B[i] * _MSX_B[j]));
                #pragma omp ordered
                {
                    B.insertBackUncompressed(j, i) =  temp;
                }
                //out_debug << j << "\t" << i << "\t" << _snp_name[_include[j]] << "\t" << _snp_name[_include[i]] << "\t" << temp << endl;
            }
        }
        if((i + 1) % 1000 == 0 || (i + 1) == _include.size()) LOGGER << i + 1 << " of " << _include.size() << " SNPs.\r";
    }
    //out_debug.close();
  
    LOGGER << "Estimating the joint effects of all SNPs ..." << endl;
    // change here
    //eigenVector Xty = D.array() * _beta.array();
    SimplicialLDLT<eigenSparseMat> solver;
    solver.compute(B);
    if(solver.info()!=Success) {
        LOGGER.e(0, "decomposition failed. The SNP correlation matrix is not positive definite.");
    }
    bJ = solver.solve(Xty);
    if(solver.info()!=Success) {
        LOGGER.e(0, "solving failed. Unable to solve the BLUP equation.");
    }

    // debug
    //LOGGER<<"_MSX: "<<_MSX<<endl;
    //LOGGER<<"Nd: "<<_Nd<<endl;
    //LOGGER<<"B: "<<B<<endl;
    //LOGGER<<"D: "<<D<<endl;
    //LOGGER<<"_beta: "<<_beta<<endl;
    //LOGGER<<"bJ: "<<bJ<<endl;

    //LOGGER<<"B inverse: "<<ldlt.solve(eigenMatrix::Identity(_include.size(), _include.size()))<<endl;


    return (true);
}



