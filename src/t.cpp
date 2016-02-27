#include "arc.h"
#include "t.h"
#include <vector>
#include <iostream>

using namespace std;

vector<Arc> gerarTransposicao(int i, int j, int k, int n){

    vector<Arc> transposicao(n);

    for(int aux = 0; aux < i; aux++){
        Arc arc(aux,aux);
        transposicao[aux] = arc;
    }

    for(int aux = i; aux < j; aux++){
        Arc arc;
        arc.a = aux;
        arc.b = aux + (k - j);
        transposicao[aux] = arc;
    }

    for(int aux = j; aux < k; aux++){
        Arc arc;
        arc.a = aux;
        arc.b = aux - (j - i);
        transposicao[aux] = arc;
    }

    for(int aux = k; aux < n; aux++){
        Arc arc(aux,aux);
        transposicao[aux] = arc;
    }

    return transposicao;
}


vector<vector<Arc>> gerarTransposicoes(int n){

    vector<vector<Arc>> transposicoes;
    int o = 0;

    vector<Arc> transposicaoNula(n);
    for(int i = 0; i < n; i++){
        Arc arc(i,i);
        transposicaoNula[i] = arc;
    }

    transposicoes.push_back(transposicaoNula);

    for(int i = 0; i <= (n-2); i++){
        for(int j = i + 1; j <= (n-1); j++){
           for(int k = j + 1; k <= n; k++){
                transposicoes.push_back(gerarTransposicao(i,j,k,n));
           }
        }
    }

    return transposicoes;
}

bool contains(Arc arc,std::vector<Arc> arcs){

    for(int i = 0; i < arcs.size(); i++){
        if(arc.Equals(arcs[i])){
            return true;
        }
    }



    return false;
}

