#pragma once

#include "gurobi_c++.h"

// Inspired by 'callback_c++.cpp' file from the GUROBI examples

class PartitionCallback : public GRBCallback {
public:
    bool print;
    int opt;
    PartitionCallback(int opt_in, bool print_in) {
        opt = opt_in;
        print = print_in;
    }
protected:
    void callback() {
        try {
            if (where == GRB_CB_POLLING) {
                // Ignore polling callback
            }
            else if (where == GRB_CB_PRESOLVE) {
                // Presolve callback

                // Ignore
            }
            else if (where == GRB_CB_SIMPLEX) {
                // Simplex callback
                /*
                double itcnt = getDoubleInfo(GRB_CB_SPX_ITRCNT);
                if (itcnt - lastiter >= 100) {
                    lastiter = itcnt;
                    double obj = getDoubleInfo(GRB_CB_SPX_OBJVAL);
                    int ispert = getIntInfo(GRB_CB_SPX_ISPERT);
                    double pinf = getDoubleInfo(GRB_CB_SPX_PRIMINF);
                    double dinf = getDoubleInfo(GRB_CB_SPX_DUALINF);
                    char ch;
                    if (ispert == 0)      ch = ' ';
                    else if (ispert == 1) ch = 'S';
                    else                  ch = 'P';
                    std::cout << itcnt << " " << obj << ch << " "
                        << pinf << " " << dinf << std::endl;
                }
                */
            }
            else if (where == GRB_CB_MIP) {
                // General MIP callback
                //double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
                //double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
                double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
                //int solcnt = getIntInfo(GRB_CB_MIP_SOLCNT);
                if (opt != -112233445566) {
                    // This condition checks whether the model has already been solved once
                    // -112233445566 is the initialization value of 'opt'
                    if (objbnd < opt) {
                        if (print) {
                            std::cout << "Stop early - objective bound lower than optimal solution" << std::endl;
                        }
                        abort();
                    }
                }
            }
            /*else if (where == GRB_CB_MIPSOL) {
                // MIP solution callback
                int nodecnt = (int)getDoubleInfo(GRB_CB_MIPSOL_NODCNT);
                double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
                int solcnt = getIntInfo(GRB_CB_MIPSOL_SOLCNT);
                double* x = getSolution(vars, numvars);
                std::cout << "**** New solution at node " << nodecnt
                    << ", obj " << obj << ", sol " << solcnt
                    << ", x[0] = " << x[0] << " ****" << std::endl;
                delete[] x;
            }
            else if (where == GRB_CB_MIPNODE) {
                // MIP node callback
                std::cout << "**** New node ****" << std::endl;
                if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {
                    double* x = getNodeRel(vars, numvars);
                    setSolution(vars, x, numvars);
                    delete[] x;
                }
            }
            else if (where == GRB_CB_BARRIER) {
                // Barrier callback
                int itcnt = getIntInfo(GRB_CB_BARRIER_ITRCNT);
                double primobj = getDoubleInfo(GRB_CB_BARRIER_PRIMOBJ);
                double dualobj = getDoubleInfo(GRB_CB_BARRIER_DUALOBJ);
                double priminf = getDoubleInfo(GRB_CB_BARRIER_PRIMINF);
                double dualinf = getDoubleInfo(GRB_CB_BARRIER_DUALINF);
                double cmpl = getDoubleInfo(GRB_CB_BARRIER_COMPL);
                std::cout << itcnt << " " << primobj << " " << dualobj << " "
                    << priminf << " " << dualinf << " " << cmpl << std::endl;
            }*/
        }
        catch (GRBException e) {
            std::cout << "Error number: " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        }
        catch (...) {
            std::cout << "Error during callback" << std::endl;
        }
    }
};
