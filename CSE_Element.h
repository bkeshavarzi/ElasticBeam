#ifndef CSE_ELEMENT_H
#define CSE_ELEMENT_H
#include <iostream>
#include "Node.h"
#include "ElasticMaterial.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>

using namespace std;
using namespace Eigen;

class CSE_Element
{
    public:

        CSE_Element();
        CSE_Element(int,double,vector <Node>,ElasticMaterial); //id,thickness,nodal objects,material
        int ldof[6];
        void SetId(int);
        int GetId(void);
        void Setth(double);
        double Getth(void);
        void SetGama(double);
        double GetGama(void);
        void SetMat(ElasticMaterial);
        ElasticMaterial GetMat(void);
        void SetNodeObj(Node,Node,Node);
        vector <Node> GetNodeObj(void);
        void SetLocalCord(void);
        MatrixXd GetLocalCord(void);
        MatrixXd SetCMatrix();
        void SetA();
        double GetA();
        void SetBMatrix();
        MatrixXd GetBMatrix(void);
        void SetDMatrix_PlaneStress();
        void SetDMatrix_PlaneStrain();
        void SetLSM();
        MatrixXd GetLSM();
        void SetU(MatrixXd);
        MatrixXd CalculateStrainTensor();
        MatrixXd CalculateStressTensor();
        MatrixXd CalculatePrinStress(string);
        MatrixXd CalculatePrinStrain(string);
        virtual ~CSE_Element();

    protected:

    private:
        int id;
        vector <Node> NodeObj;
        ElasticMaterial mat;
        double E,v,gama,th,A;
        MatrixXd LCord=MatrixXd::Zero(2,3);
        MatrixXd LSM=MatrixXd::Zero(6,6);
        MatrixXd D=MatrixXd::Zero(3,3);
        MatrixXd U=MatrixXd::Zero(6,1);
        MatrixXd Ep=MatrixXd::Zero(3,1);
        MatrixXd Sigma=MatrixXd::Zero(3,1);
        MatrixXd PSigma=MatrixXd::Zero(3,1);
        MatrixXd PStrain=MatrixXd::Zero(3,1);
};

#endif // CSE_ELEMENT_H
