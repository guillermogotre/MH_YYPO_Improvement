//
// Created by guillermo on 6/22/18.
//

#include "YYOPMemeDifStop.h"
#include "localsearch.h"
#include "cmaeshan.h"
#include "solis.h"
#include "simplex.h"

#include <algorithm>
#include <math.h>
#include <assert.h>
#include <float.h>

tFitness realea::YYOPMemeDifStop::evalSol(tChromosomeReal &orig) {
    tChromosomeReal fitted(orig.size());
    std::transform(orig.begin(),orig.end(),fitted.begin(),
                   [&](tReal v){return v*(maxV-minV)+minV;}
    );
    return m_problem->eval(fitted);
}
realea::tChromosomeReal realea::YYOPMemeDifStop::solToDim(tChromosomeReal & s){
    tChromosomeReal a(s.size());
    std::transform(s.begin(),s.end(),a.begin(),
                   [&](tReal v){return v*(maxV-minV)+minV;}
    );
    return a;
}

realea::tChromosomeReal realea::YYOPMemeDifStop::solToOne(tChromosomeReal & s){
    tChromosomeReal a(s.size());
    std::transform(s.begin(),s.end(),a.begin(),
                   [&](tReal v){return (v-minV)/(maxV-minV);}
    );
    return a;
}

tFitness realea::YYOPMemeDifStop::splitStage(tChromosomeReal & s, tReal d){
    const double SPLIT_PROB = 0.5;
    const unsigned D = s.size();
    const double DWAYFACTOR = sqrt(2);

    vector<tChromosomeReal> S(2*D,s);

    double R = m_random->randreal(0,1);
    //One-way splitting
    if(R > SPLIT_PROB){
        for(unsigned i=0;i<s.size();i++){
            //I
            tReal newR  = m_random->randreal(0,1)*d + s[i];
            if(newR > 1 || newR < 0)
                newR = m_random->randreal(0,1);
            S[i][i] = newR;
            //D+I
            newR  = m_random->randreal(0,1)*d - s[i];
            if(newR > 1 || newR < 0)
                newR = m_random->randreal(0,1);
            S[D+i][i] = newR;
        }
        //D-way splitting
        //TODO not equal
    }else{
        for(unsigned i=0; i<S.size(); i++){
            for(unsigned j=0; j<D; j++){
                tReal newR;
                if(m_random->randint(0,1)==1){
                    newR  = m_random->randreal(0,1)*(d/DWAYFACTOR) + s[j];
                }else{
                    newR  = m_random->randreal(0,1)*(d/DWAYFACTOR) - s[j];
                }
                if(newR > 1 || newR < 0)
                    newR = m_random->randreal(0,1);
                S[i][j] = newR;
            }
        }
    }
    int bestIndex = -1;
    tFitness bestFitness = DBL_MAX;
    for(unsigned i=0; i<S.size(); i++){
        tFitness f = evalSol(S[i]);
        if(f<bestFitness){
            bestFitness = f;
            bestIndex = i;
        }
    }
    assert(bestIndex != -1);
    s = S[bestIndex];
    return bestFitness;
}

unsigned realea::YYOPMemeDifStop::lsStage(tChromosomeReal & sol, tFitness & fitness, unsigned iters){
    /*
     * Local Search
     * */
    ILocalSearch *ls;
    ILSParameters *ls_options;

    ls = new SimplexDim();

    ls->setProblem(m_problem);
    ls->setRandom(m_random);

    tChromosomeReal p = solToDim(sol);

    ls_options = ls->getInitOptions(p);

    return ls->apply(ls_options,p,fitness,iters);
}

unsigned realea::YYOPMemeDifStop::difStage(vector<tChromosomeReal> &arc, vector<tFitness> & fitness) {
    const unsigned POP_SIZE = arc.size();
    const unsigned SOL_SIZE = arc[0].size();
    const double CR = 0.5, F=0.5, LS = 0.1;
    const double LSFACTOR = 0.5;

    unsigned evals = 0;

    vector<tChromosomeReal> mut = arc;
    vector<tFitness> evs(POP_SIZE);

    for(int i=0; i<POP_SIZE; i++){
        //Select 3 random parents
        //cout << ";" << endl;
        int     r1 = m_random->randint(0,(POP_SIZE-1)/2)*2,
                r2 = m_random->randint(0,(POP_SIZE-1)/2)*2,
                r3 = m_random->randint(0,(POP_SIZE-1)/2)*2 + 1;
        while(r2 == r1 || r2 == i) r2 = m_random->randint(0,POP_SIZE-1);
        while(r3 == r2 || r3 == r1 || r3 == i) r3 = m_random->randint(0,POP_SIZE-1);
        //Mutate
        //For gene
        for(int j=0; j<SOL_SIZE; j++){
            //cout << "1" << endl;
            double r = m_random->randreal(0,1);
            if(r<CR){
                //cout << "2:" << r1 << " " << r2 << " " << r3 << endl;
                mut[i][j] = arc[r1][j] + F*(arc[r2][j]-arc[r3][j]);
                //cout << "3" << endl;
                //Check bound
                if (mut[i][j] > 1)
                    mut[i][j] = 1;
                else if(mut[i][j] < 0)
                    mut[i][j] = 0;
            }
        }
        evs[i] = evalSol(mut[i]);
        double r = m_random->randreal(0,1);
        if(r<LS){
            //cout << "#" << evs[i] << endl;
            evals += lsStage(mut[i],evs[i],LSFACTOR*SOL_SIZE);
            //cout << evs[i] << endl;
        }
        evals++;
    }
    for(int i=0; i<POP_SIZE; i++){
        if(evs[i] < fitness[i]){
            arc[i] = mut[i];
            fitness[i] = evs[i];
        }
    }
    return evals;
}



unsigned realea::YYOPMemeDifStop::apply(tChromosomeReal &sol, tFitness &fitness, unsigned itera) {

    const unsigned I_MIN = 2, I_MAX = 10, I = m_random->randint(I_MIN,I_MAX);
    const unsigned MAX_ITER = itera;
    const double MAX_D = 0.75;


    auto getInitRandom = [&](tReal min, tReal max){
        const  int D = m_problem->getDimension();
        tChromosomeReal p(D);
        for(unsigned i=0; i<D; i++){
            p[i] = m_random->randreal(min,max);
        }
        return p;
    };

    //cout << "I: " << I << endl;
    int iters = itera;
    tReal alpha = 10;
    tReal d1 = 0.5, d2 = 0.5;
    tChromosomeReal p1 = solToOne(sol), p2 = getInitRandom(0,1);



    vector<tChromosomeReal> archive(2*I);
    vector<tReal> archiveEval(2*I);

    tChromosomeReal fittedSol(sol.size());

    tReal eval1 = evalSol(p1);
    tReal eval2 = evalSol(p2);

    /*
     * STOP MUTATION
     * */
    const int NOBETTER = p1.size()*3;
    const double MUT_PROB = 0.8;
    int failCount = 0;
    tChromosomeReal bestSol = p1;
    tReal bestHistEval = eval1;
    tReal bestCurrEval = eval1;
    tChromosomeReal initSol = p1;
    tFitness  initEval = eval1;

    int count1 = 0, count2 = 0;
    iters-=2;

    tReal bestEval = DBL_MAX;

    //Aux functions
    auto swp = [&](){p1.swap(p2);swap(eval1,eval2),swap(d1,d2);swap(count1,count2);};


    int i=0;
    int LS_ITERS = 0.1*MAX_ITER;
    int LS_ITERS1 = 0.6*LS_ITERS;
    int LS_ITERS2 = 1.4*MAX_ITER;
    int added = p1.size()*4*I;

    /*
     * Local Search
     * */
    ILocalSearch *ls;
    ILSParameters *ls_options;

    CMAESHansen *cmaes = new CMAESHansen("cmaesinit.par");
    cmaes->searchRange(0.1);
    ls = cmaes;


    ls->setProblem(m_problem);
    ls->setRandom(m_random);

    ls_options = ls->getInitOptions(p1);
    p1 = solToDim(p1);

    iters -= ls->apply(ls_options,p1,eval1,LS_ITERS1);
    p1 = solToOne(p1);

    vector<unsigned> indices(archive.size());
    unsigned n = 0;
    generate(indices.begin(),indices.end(),[&]{ return n++;});


    iters += MAX_ITER/4;
    while((iters-added-LS_ITERS*2) > 0 || i != 0){ //After archive stage
        //while(false){
        //Check fittest
        if(eval2 < eval1){
            swp();
            //cout << MAX_ITER-iters << " " << i << " " << d1 << " " << d2 << " SWAP!" << " ["<< count1 <<"," <<  count2<<"] "<< endl;
        }
        if(eval1 < bestEval)
            bestEval = eval1;

        //Add to archive
        unsigned offset = i*2;

        archive[offset] = p1;
        archive[offset+1] = p2;
        archiveEval[offset] = eval1;
        archiveEval[offset+1] = eval2;

        i++;
        count1++;

        //Archive stage
        if(i == I){
            i=0;
            iters -= difStage(archive,archiveEval);
            //Select two fittest
            sort(indices.begin(),indices.end(),
                 [&](int i1, int i2){
                     return archiveEval[i1] < archiveEval[i2];
                 });
            p1 = archive[indices[0]];
            eval1 = archiveEval[indices[0]];

            p2 = archive[indices[1]];
            eval2 = archiveEval[indices[1]];

            //Update d
            d1 = d1 - (d1/alpha);
            d2 = min(d2 + (d2/alpha),MAX_D);


            //if no better
            if(eval1<bestCurrEval){
                bestCurrEval = eval1;
                if(eval1 < bestHistEval){
                    bestHistEval = eval1;
                    bestSol = p1;
                }
            }
            else{
                failCount++;
                if(failCount > NOBETTER){
                    //cout << "RESET" << endl;
                    failCount = 0;
                    d1 = 0.5;
                    d2 = 0.5;
                    //Mutate p1
                    for(int j=0; j<p1.size(); j++){
                        if(m_random->randreal(0,1)<0.5)
                            p1[j] = p1[j]+m_random->normal(0.3);
                        if(p1[j] < 0 | p1[j] > 1)
                            p1[j] = m_random->randreal(0,1);
                    }
                    eval1 = evalSol(p1);
                    iters -= lsStage(p1,eval1,p1.size()*5);
                    //Random p2
                    p2 = getInitRandom(0,1);
                    eval2 = evalSol(p2);
                    bestCurrEval = eval1;

                    iters -= 2;
                }
            }
        }
            //Splitting stage
        else{
            eval1 = splitStage(p1,d1);
            eval2 = splitStage(p2,d2);

            iters -= 4*p1.size();
        }
    }
    //assert(bestEval == eval1);
    //printDif(p1,p2);
    p1 = solToDim(bestSol);
    cmaes->searchRange(0.0004);
    //cmaes->searchNeighborhood(0.01,&archive);
    ls_options = ls->getInitOptions(p1);
    iters -= ls->apply(ls_options,p1,bestHistEval,LS_ITERS2);

    tFitness before,after,diff;

    before = fitness;

    //cout << "Evals: " << MAX_ITER - iters << endl;

    after = bestHistEval;

    diff = before-after;

    //cout <<"Improvement: " <<std::scientific <<before <<" -> " <<std::scientific <<after;
    //cout <<" [" <<std::scientific <<diff <<"]" << " , Optimo: [" << m_problem->getOptime() << "]" << endl;

    sol = bestSol;
    fitness = bestHistEval;
}