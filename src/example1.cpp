#include <iostream>
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

/*int main(int argc, char** argv)
{
    std::cout << "Iniciando exemplo 1" << endl;

    IloEnv env;

    try{

        IloModel model(env);

        IloNumVarArray vars(env);
        vars.add(IloNumVar(env,0.0,40.0));
        vars.add(IloNumVar(env));
        vars.add(IloNumVar(env));
        model.add(IloMaximize(env, vars[0] + 2 * vars[1] + 3 * vars[2]));
        model.add(-vars[0] + vars[1] + vars[2] <= 20);
        model.add(vars[0] - vars[1] + vars[2] <= 30);

        IloCplex cplex(model);

        if(!cplex.solve()){
            env.error() << "Failed to optimize LP" <<endl;
            throw(-1);
        }

        IloNumArray vals(env);
        env.out() << "Solution status = " << cplex.getStatus() << endl;
        env.out() << "Solution valur = " << cplex.getObjValue() << endl;
        cplex.getValues(vals,vars);
        env.out() << "Values = " << vals << endl;

    }catch(IloException& e){
        cerr << "Concert exception caught: " << e << endl;
    }catch(...){
        cerr << "Unknow exception caught" << endl;
    }

    env.end();

    return 0;
}*/
