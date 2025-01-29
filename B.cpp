#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "glpk.h"

using namespace std;

int n, m;
vector<pair<int, int> > adjVerices;
vector<vector<int> > adjEdges;
int comp[500];
int curCompSize;

glp_prob *mip;

void dfs(const int v, const int c_num, const vector<vector<int> > &adj) {
    comp[v] = c_num;
    curCompSize++;

    for (const int u: adj[v]) {
        if (comp[u] == 0) {
            dfs(u, c_num, adj);
        }
    }
}

void add_constraints(const int oddCnt, const unordered_map<int, int> &oddComp) {
    const int fstRow = glp_add_rows(mip, oddCnt);
    unordered_map<int, vector<int> > eMap;
    for (const auto [cNum, cSize]: oddComp) {
        eMap.insert({cNum, vector<int>()});
    }
    for (int e = 0; e < m; e++) {
        if (const int cNum = comp[adjVerices[e].first];
            cNum == comp[adjVerices[e].second] && oddComp.find(cNum) != oddComp.end()) {
            eMap[cNum].push_back(e + 1);
        }
    }
    int rIdx = 0;
    for (const auto [cNum, cSize]: oddComp) {
        vector<int> es = eMap[cNum];
        const int ne = es.size();
        double val[ne + 1];
        int ids[ne + 1];
        fill_n(val, ne + 1, 1.0);
        std::copy(es.begin(), es.end(), ids + 1);

        glp_set_row_bnds(mip, fstRow + rIdx, GLP_UP, 0.0, 0.5 * (cSize - 1));
        glp_set_mat_row(mip, fstRow + rIdx, ne, ids, val);
        rIdx++;
    }
}

void gen_false_constraints(glp_tree *T) {
    auto adj = vector(n, vector<int>());
    fill_n(comp, n, 0);
    for (int e = 0; e < m; e++) {
        if (glp_get_col_prim(mip, e + 1) >= .001) {
            adj[adjVerices[e].first].push_back(adjVerices[e].second);
            adj[adjVerices[e].second].push_back(adjVerices[e].first);
        }
    }

    int cNum = 1;
    unordered_map<int, int> oddComp;
    int oddCnt = 0;
    for (int i = 0; i < n; i++) {
        if (comp[i] == 0) {
            curCompSize = 0;
            dfs(i, cNum, adj);
            if (curCompSize % 2 == 1 && curCompSize > 1) {
                oddCnt++;
                oddComp.insert({cNum, curCompSize});
            }
            cNum++;
        }
    }

    if (oddCnt > 0) {
        add_constraints(oddCnt, oddComp);
    }
}

void cb_func(glp_tree *T, void *info) {
    switch (glp_ios_reason(T)) {
        case GLP_IROWGEN:
            gen_false_constraints(T);
            break;
    }
}

int main() {
     ifstream cin("min-weight-perfect-matching/01");

    ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    cin >> n >> m;

    double edgeCost[m];
    adjEdges = vector(n, vector<int>());
    for (int e = 0; e < m; e++) {
        int a, b;
        cin >> a >> b >> edgeCost[e];
        adjEdges[a].push_back(e);
        adjEdges[b].push_back(e);
        adjVerices.emplace_back(a, b);
    }
    int matrixN = m * 2;

    int ia[matrixN + 1], ja[matrixN + 1];
    double ar[matrixN + 1];
    int idxA = 1;

    mip = glp_create_prob();
    glp_set_prob_name(mip, "matching");
    glp_set_obj_dir(mip, GLP_MIN);

    glp_add_rows(mip, n);
    for (int i = 0; i < n; i++) {
        glp_set_row_name(mip, i + 1, ("v" + to_string(i)).c_str());
        glp_set_row_bnds(mip, i + 1, GLP_FX, 1.0, 1.0);
    }

    glp_add_cols(mip, m);
    for (int e = 0; e < m; e++) {
        glp_set_col_name(mip, e + 1, ("x" + to_string(e)).c_str());
        glp_set_col_kind(mip, e + 1, GLP_BV);
        glp_set_obj_coef(mip, e + 1, edgeCost[e]);
    }

    for (int v = 0; v < n; v++) {
        for (const int e: adjEdges[v]) {
            ia[idxA] = v + 1;
            ja[idxA] = e + 1;
            ar[idxA] = 1.0;
            idxA++;
        }
    }

    glp_load_matrix(mip, idxA - 1, ia, ja, ar);

    glp_smcp parm1;
    glp_init_smcp(&parm1);
    parm1.msg_lev = GLP_MSG_OFF;
    glp_simplex(mip, &parm1);

    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.cb_func = cb_func;
    parm.msg_lev = GLP_MSG_OFF;
    parm.gmi_cuts = GLP_ON;
    glp_intopt(mip, &parm);

    const double z = glp_mip_obj_val(mip);
    cout << z << endl;

    for (int e = 0; e < m; e++) {
        if (const double y = glp_mip_col_val(mip, e + 1); y > 0.5) {
            cout << e << " ";
        }
    }
    cout << endl;

    glp_delete_prob(mip);
    return 0;
}
