#ifndef TEXTFILES_H_INCLUDED
#define TEXTFILES_H_INCLUDED
#include <iostream>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include "Node.h"
#include "CSE_Element.h"
#include "Q4_Element.h"
#include "Q9_Element.h"
#include "Q8_Element.h"

using namespace std;
using namespace Eigen;

vector <Node> ReadNodeFile(string);
vector <CSE_Element> ReadCSEElement(string,ElasticMaterial &,double,vector <Node>);
vector <Q4_Element>  ReadQ4Element(string,ElasticMaterial &,double,vector <Node>);
vector <Q8_Element>  ReadQ8Element(string,ElasticMaterial &,double,vector <Node>);
vector <Q9_Element>  ReadQ9Element(string,ElasticMaterial &,double,vector <Node>);
vector <Node> SortElementNodeQ4(vector <Node>);
vector <Node> SortElementNodeQ8(vector <Node>);
vector <Node> SortElementNodeQ9(vector <Node>);
void WriteOutPutFile_CSE(MatrixXd, vector <Node>,vector <CSE_Element>);
//MatrixXd MakeNodeMatrix(vector <Node>);
//MatrixXd MakeQ4Matrix(vector <Q4_Element>);
//MatrixXd MakeQ8Matrix(vector <Q8_Element>);
//MatrixXd MakeQ9Matrix(vector <Q9_Element>);

#endif // TEXTFILES_H_INCLUDED
