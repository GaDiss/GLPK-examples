#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "glpk.h"

using namespace std;

int main() {
     ifstream cin("capacitated-warehouse-location/01");

    ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n, m;
    cin >> n >> m;

    double capacity[n+1];
    double openCost[n+1];
    for (int w = 1; w <= n; w++) {
        cin >> capacity[w] >> openCost[w];
    }

    double demand[m+1];
    for (int c = 1; c <= m; c++) {
        cin >> demand[c];
    }

    double useCost[n+1][m+1];
    for (int w = 1; w <= n; w++) {
        for (int c = 1; c <= m; c++) {
            cin >> useCost[w][c];
        }
    }

    int rowsN = m + n;
    int colsN = n + n * m;
    int matrixN = n + (n * m) * 2;

    int ia[matrixN+1], ja[matrixN+1];
    double ar[matrixN+1];
    int idxA = 1;

    glp_prob *mip = glp_create_prob();
    glp_set_prob_name(mip, "warehouses");
    glp_set_obj_dir(mip, GLP_MIN);

    glp_add_rows(mip, rowsN);
    for (int c = 1; c <= m; c++) {
        glp_set_row_name(mip, c, ("cr" + to_string(c)).c_str());
        glp_set_row_bnds(mip, c, GLP_FX, 1.0, 1.0);
    }
    for (int w = 1; w <= n; w++) {
        glp_set_row_name(mip, m + w, ("wr" + to_string(w)).c_str());
        glp_set_row_bnds(mip, m + w, GLP_UP, 0.0, 0.0);
    }

    glp_add_cols(mip, colsN);
    for (int w = 1; w <= n; w++) {
        glp_set_col_name(mip, w, ("y" + to_string(w)).c_str());
        glp_set_col_kind(mip, w, GLP_BV);
        glp_set_obj_coef(mip, w, openCost[w]);

        ia[idxA] = m+w;
        ja[idxA] = w;
        ar[idxA] = -capacity[w];
        idxA++;
    }

    for (int c = 1; c <= m; c++) {
        for (int w = 1; w <= n; w++) {
            glp_set_col_name(mip, c * n + w, ("x" + to_string(w) + "." + to_string(c)).c_str());
            glp_set_col_bnds(mip, c * n + w, GLP_LO, 0.0, 0.0);
            glp_set_obj_coef(mip, c * n + w, useCost[w][c]);

            ia[idxA] = c;
            ja[idxA] = c * n + w;
            ar[idxA] = 1.0;
            idxA++;

            ia[idxA] = m+w;
            ja[idxA] = c * n + w;
            ar[idxA] = demand[c];
            idxA++;

        }
    }

     glp_load_matrix(mip, idxA-1, ia, ja, ar);

    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    parm.msg_lev = GLP_MSG_OFF;
    // parm.tm_lim = 19900;
    parm.mip_gap = 0.001;
    parm.gmi_cuts = GLP_ON;
    parm.mir_cuts = GLP_ON;
    glp_intopt(mip, &parm);

    int openW = 0;
    vector<int> ow;
    for (int w= 1; w <= n; w++) {
        if (double y = glp_mip_col_val(mip, w); y > 0.5) {
            openW+=1;
            ow.push_back(w);
        }
    }
    cout << openW << endl;
    for (int w : ow) {
        cout << w << " ";
    }
    cout << endl;


    cout << fixed << showpoint;
    for (int w : ow) {
        for (int c = 1; c <= m; c++) {
            double x = glp_mip_col_val(mip, c * n + w);
            cout << std::setprecision(6) << x << ' ';
        }
        cout << endl;
    }

    // double z = glp_mip_obj_val(mip);
    // cout << std::setprecision(6) << z << endl;

    glp_delete_prob(mip);
    return 0;
}
