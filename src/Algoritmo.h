//
// Created by guillermo on 6/11/18.
//

#ifndef LS_ALGORITMO_H
#define LS_ALGORITMO_H

#include "problem.h"
#include "random.h"

namespace realea {

    class Algoritmo {
    protected:
        Random *m_random;
        IEval *m_eval;
        Problem *m_problem;
        //Every problem is defined in the same
        //limit for every dimension (-+100)
        tReal minV, maxV;

        tFitness evalSol(tChromosomeReal & orig);
        tFitness splitStage(tChromosomeReal &s, tReal d);

    public:
        unsigned apply(tChromosomeReal & sol, tFitness & fitness, unsigned itera);

        void setProblem(Problem * problem){
            m_problem = problem;
            m_problem->getDomain()->getValues(0,&minV,&maxV);
            setEval(problem);
        }

        void setRandom(Random * random){
            m_random = random;
        }

        void setEval(IEval * eval){
            m_eval = eval;
        }

    };
}
#endif //LS_ALGORITMO_H
