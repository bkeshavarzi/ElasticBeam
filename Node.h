#ifndef NODE_H
#define NODE_H
#include <iostream>
#include <vector>

using namespace std;

class Node
{
    public:

        Node();
        Node(int,double,double,double,double,bool,bool);
        void SetId(int);
        int GetId();
        void SetCord(double,double);
        vector <double> GetCord(void);
        void SetForce(double,double);
        vector <double> GetForce(void);
        void SetRes(bool,bool);
        vector <bool> GetRes(void);
        virtual ~Node();

    protected:

    private:

        int id;
        double xc,yc;
        double fx,fy;
        bool rx,ry;
};

#endif // NODE_H
