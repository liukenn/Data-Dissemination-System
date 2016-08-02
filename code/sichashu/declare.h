#include<iostream>
#include<stdlib.h>
#include<cstdlib>
#include<vector>
#include<fstream>
#include<cmath>
#include<ctime>
#include <algorithm>
using namespace std;



//空间对象为平面直角坐标系里一个坐标
//象限说明
//    2    |    3
//         |    
//---------|-----------
//         |
//    1    |    0

typedef struct ObInf				//目标对象坐标结构体
{
	int l;
	int r;
}ObInf;

typedef struct QuadBox				//表示范围的结构体
{
	int maxl;
	int minl;
	int maxr;
	int minr;
}QuadBox;

typedef struct QuadNode				//四叉树节点结构体
{
	QuadBox box;
	vector<ObInf> object;
	int ObCount;
	struct QuadNode *parent;
	struct QuadNode *children[4];
}QuadNode;

typedef struct QuadTree				//四叉树结构体
{
	QuadNode *p;
	int depth;
}QuadTree;

typedef struct SplitPointInf
{
	double d;
	ObInf coordinate;
}SplitPointInf;


QuadTree *InitQuadTree(QuadBox box,vector<vector<int> > AdjMatrix);
void AddObject(ObInf splitpoint,QuadNode *p,ObInf ob);
QuadNode *AddNode(QuadNode *p,QuadBox box,QuadNode *parent);
QuadBox &InitQuadBox(vector<vector<int> > AdjMatrix);
void ViewTree(QuadNode *p,int n);
int GetNodeHight(QuadNode *p);
int &CalculateTreeHight(const int v);
vector<vector<int> > CreatDenMatrix(vector<vector<int> > AdjMatrix);
ObInf &GetSplitPoint(vector<vector<int> > AdjMatrix,vector<vector<int> > DenMatrix,int minsize,int splittimes,int nh,int v);
void SplitQuadNode(QuadTree *pt,QuadNode *p,const int splittimes,vector<vector<int> > AdjMatrix);
void InitQuadNode(ObInf splitpoint,QuadNode *p);
int NoisyCount(//double epsilon,				//隐私预算
			   double ranNum,
			   int ObCount,					//叶子节点中1的个数
			   int i,						//当前树的高度
			   int h						//树的高度
	);
double RanNumber(int n);
vector<double> BubbleSort(vector<double> L,int n);//冒泡排序
int ExpCount(int i,//当前树的高度
			 int h,//树的高度
			 int v,//矩阵大小
			 int n,//分裂点个数，分裂点数组的大小
			 vector<double> dis//存储分裂点对应密度差
			 );