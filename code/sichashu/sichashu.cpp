
#include"stdafx.h"


#include"declare.h"
using namespace std;

int main()
{
	vector<vector<int> > AdjMatrix,DenMatrix;									//AdjMatrix：存储邻接矩阵 DenMatrix：存储生成的密度矩阵 都使用二维的vector存储

	int a;
	ifstream fin("adjmatrixtest.txt");
	//ifstream fin("D://test(4).txt");
	if(fin==NULL)
	{
		cout<<"open file failed!"<<endl;
		return 0;
	}
	int v;
	cout<<"输入邻接矩阵大小v："<<endl;											//矩阵大小需要手动输入
	cin>>v;
	vector<int> temp(v);
	for(int i=0;i<v;++i)														//读取存储在测试文件中的邻接矩阵 实际程序中之前应该已经建立了 所以没写成函数
	{
		for(int j=0;j<v;++j)
		{
			fin>>a;
			temp[j]=a;
		}
		AdjMatrix.push_back(temp);

	}																			//建立邻接矩阵
	cout<<"已成功创建邻接矩阵与密度矩阵！"<<endl;								//接下来使用这个邻接矩阵建立四叉树


	QuadTree *pt;																//建立与分裂四叉树
	QuadBox box;
	vector<ObInf> ob;
	box=InitQuadBox(AdjMatrix);													
	pt=InitQuadTree(box,AdjMatrix);
	cout<<"已根据邻接矩阵初始化四叉树!"<<endl;
	int splittimes=CalculateTreeHight(v);	
	SplitQuadNode(pt,pt->p,splittimes,AdjMatrix);
	cout<<"分裂完的四叉树高度是"<<pt->depth<<endl;
	int n=1;
	ViewTree(pt->p,n);	//输出分裂完的四叉树

	return 0;
}



QuadBox &InitQuadBox(vector<vector<int> > AdjMatrix)
{
	int size;
	QuadBox box;
	size=AdjMatrix.size();
	box.maxl=size;
	box.maxr=size;
	box.minl=1;
	box.minr=1;
	return box;
}


QuadTree *InitQuadTree(QuadBox box,vector<vector<int> > AdjMatrix)					//初始化四叉树，需要参数为范围大小box与之前生成的邻接矩阵
{
	QuadTree *pt;
	pt=(QuadTree *)new(QuadTree);
	if(pt==NULL)
	{
		cout<<"初始空间分配错误"<<endl;
		system("PAUSE");
		return 0;		
	}
	pt->p=AddNode(pt->p,box,NULL);
	pt->depth=0;
	int size=AdjMatrix.size();
	ObInf tempob;
	for(int i=0;i<size;++i)
		for(int j=0;j<size;++j)
		{
			if(AdjMatrix[i][j]==1)
			{
				tempob.l=i+1;
				tempob.r=j+1;
				pt->p->object.push_back(tempob);
				++(pt->p->ObCount);
			}
		}
	return pt;
}


QuadNode *AddNode(QuadNode *p,QuadBox box,QuadNode *parent)							//添加节点
{
	p=(QuadNode *)new(QuadNode);
	if(p!=NULL)
	{
		p->ObCount=0;
		p->box.maxl=box.maxl;
		p->box.maxr=box.maxr;
		p->box.minl=box.minl;
		p->box.minr=box.minr;
		p->children[0]=NULL;
		p->children[1]=NULL;
		p->children[2]=NULL;
		p->children[3]=NULL;
		p->parent=parent;
		return p;
	}
	else 
	{
		cout<<"空间分配错误"<<endl;
		system("PAUSE");
		return 0;
	}
}



void AddObject(ObInf splitpoint,QuadNode *p,ObInf ob)
{
	if(ob.l>splitpoint.l&&ob.r>splitpoint.r)
	{
		p->children[0]->object.push_back(ob);
		++p->children[0]->ObCount;
	}
	else  if(ob.l>splitpoint.l&&ob.r<=splitpoint.r)
	{
		p->children[1]->object.push_back(ob);
		++p->children[1]->ObCount;
	}
	else  if(ob.l<=splitpoint.l&&ob.r<=splitpoint.r)
	{
		p->children[2]->object.push_back(ob);
		++p->children[2]->ObCount;
	}
	else  if(ob.l<=splitpoint.l&&ob.r>splitpoint.r)
	{
		p->children[3]->object.push_back(ob);
		++p->children[3]->ObCount;
	}
	return;
}



void ViewTree(QuadNode *p,int n)													//显示输出信息
{
	QuadNode *vp;
	vp=p;

	if(vp->children[0]==NULL&&vp->children[1]==NULL&&vp->children[2]==NULL&&vp->children[3]==NULL)
		cout<<"["<<vp->box.minl<<","<<vp->box.maxl<<";"<<vp->box.minr<<","<<vp->box.maxr<<"]"<<"	"<<vp->ObCount<<"--"<<
		NoisyCount(RanNumber(n),vp->ObCount,GetNodeHight(p),CalculateTreeHight(100))<<endl;
	for(int i=0;i<4;++i)
		if(vp->children[i]!=NULL)
		{
			n++;
			ViewTree(vp->children[i],n);
			
		}
	
	return;
}




int GetNodeHight(QuadNode *p)												//获得当前节点高度
{
	int h=0;
	while(p->parent!=NULL)
	{
		++h;
		p=p->parent;
	}
	return h;
}									



int &CalculateTreeHight(const int v)									//计算树的高度
{
	int h=1;
	double t=(pow(2,1/3.0)-1)*v*v/pow(2,0.5);
	while(pow(2,2*h+1/3.0)-pow(2,5*h/3.0)<=t)
		++h;
	h--;
	return h;
}



vector<vector<int> > CreatDenMatrix(vector<vector<int> > AdjMatrix)		//生成密度矩阵的函数 参数为已生成的邻接矩阵 返回生成的密度矩阵
{
	int l=AdjMatrix.size(),r=AdjMatrix[0].size();
	vector<vector<int> > DenMatrix(l,r);
	for(int i=0;i<l;++i)
	{
		for(int j=0;j<r;++j)
		{
			if(i==0&&j==0)
				DenMatrix[i][j]=AdjMatrix[i][j];
			else if(i==0&&j!=0)
				DenMatrix[i][j]=DenMatrix[i][j-1]+AdjMatrix[i][j];
			else if(i!=0&&j==0)
				DenMatrix[i][j]=DenMatrix[i-1][j]+AdjMatrix[i][j];
			else
				DenMatrix[i][j]=DenMatrix[i-1][j]+DenMatrix[i][j-1]-DenMatrix[i-1][j-1]+AdjMatrix[i][j];
		}
	}
	return DenMatrix;
}



ObInf &GetSplitPoint(vector<vector<int> > AdjMatrix,vector<vector<int> > DenMatrix,int minsize,int splittimes,int nh,int v)							//计算分裂点函数，参数为需要计算分裂点的邻接矩阵和分裂的最小区域大小，返回为（子矩阵的）分裂点坐标
{
	//int h=splittimes;//树的高度
	
	
	
	int l=AdjMatrix.size(),r=AdjMatrix[0].size();
	SplitPointInf p;
	p.d=-1;
	p.coordinate.l=-1;
	p.coordinate.r=-1;

	//double *dis=new double[]; 
	vector<double> dis;
	//int k=0;//dis[]的大小
	vector<int> larray;
	vector<int> rarray;




	for(int i=0;i<l;++i)
	{
		for(int j=0;j<r;++j)
		{
			if((i+1)*(j+1)>minsize&&(i+1)*(r-j-1)>minsize&&(l-i-1)*(j+1)>minsize&&(l-i-1)*(r-j-1)>minsize)
			{
				double d0=(DenMatrix[i][j]+0.0)/((i+1.0)*(j+1.0)),
					d1=(DenMatrix[i][r-1]-DenMatrix[i][j]+0.0)/((r-j-1.0)*(i+1.0)),
					d2=(DenMatrix[l-1][j]-DenMatrix[i][j]+0.0)/((j+1.0)*(l-i-1.0)),
					d3=(DenMatrix[l-1][r-1]-DenMatrix[l-1][j]-DenMatrix[i][r-1]+DenMatrix[i][j]+0.0)/((r-j-1.0)*(l-i-1.0));
				double d=max(d0,max(d1,max(d2,d3)))-min(d0,min(d1,min(d2,d3)));

				dis.push_back(d);
				larray.push_back(i+1);
				rarray.push_back(j+1);
				if(i>=(l/2)&&j>=(r/2)&&(l-i-1)*(r-j-2)==minsize)
				{
					int point=ExpCount(nh,splittimes,v,dis.size(),dis);
					p.d=dis[point];
					p.coordinate.l=larray[point];
					p.coordinate.r=rarray[point];
				}

				
			}
		}
	}

			
	return p.coordinate;
}	





void SplitQuadNode(QuadTree *pt,QuadNode *p,const int splittimes,vector<vector<int> > AdjMatrix)								//分裂四叉树函数
{
	int h=GetNodeHight(p);
	if(h<splittimes)
	{
		int v=AdjMatrix.size();
		int minsize=v*v/(int)pow(4.0,h+2);
		vector<vector<int> > newmatrix;
		vector<int> temp(p->box.maxr-p->box.minr+1);
		for(int i=p->box.minl-1;i<p->box.maxl;++i)
		{
			int t=0;
			for(int j=p->box.minr-1;j<p->box.maxr;++j)
			{
				temp[t]=AdjMatrix[i][j];
				++t;
			}
			newmatrix.push_back(temp);
		}
		vector<vector<int> > DenMatrix=CreatDenMatrix(newmatrix);
		ObInf point=GetSplitPoint(newmatrix,DenMatrix,minsize,splittimes,h,v);
		if(point.l==-1)
		{
			for(int i=0;i<4;++i)
				if(p->ObCount<0.8*v*v/pow(4.0,splittimes))
				{
					free(p);
					p->children[i]=NULL;
				}
				return;
		}
		ObInf npoint;
		npoint.l=point.l+p->box.minl-1;
		npoint.r=point.r+p->box.minr-1;
		InitQuadNode(npoint,p);
		for(int i=p->box.minl-1;i<p->box.maxl;++i)
			for(int j=p->box.minr-1;j<p->box.maxr;++j)
				if(AdjMatrix[i][j]==1)
				{
					ObInf ob;
					ob.l=i+1;
					ob.r=j+1;
					AddObject(npoint,p,ob);
				}
		for(int i=0;i<4;++i)
		{
			QuadNode *np;
			np=p->children[i];
			h=GetNodeHight(np);
			if(h>pt->depth)
				pt->depth=h;
			double d;
			QuadBox nbox;
			nbox.maxl=np->box.maxl-p->box.minl+1;
			nbox.minl=np->box.minl-p->box.minl+1;
			nbox.maxr=np->box.maxr-p->box.minr+1;
			nbox.minr=np->box.minr-p->box.minr+1;
			if(nbox.minl==1&&nbox.minr!=1)
				d=(DenMatrix[nbox.maxl-1][nbox.minr-1]-DenMatrix[nbox.maxl-1][nbox.minr-2])/((np->box.maxl-np->box.minl+1)*(np->box.maxr-np->box.minr+1.0));
			else if(nbox.minr==1&&nbox.minl!=1)
				d=(DenMatrix[nbox.maxl-1][nbox.minr-1]-DenMatrix[nbox.minl-2][nbox.maxr-1])/((np->box.maxl-np->box.minl+1)*(np->box.maxr-np->box.minr+1.0));
			else if(nbox.minr==1&&nbox.minl==1)
				d=DenMatrix[nbox.maxl-1][nbox.minr-1]/((np->box.maxl-np->box.minl+1)*(np->box.maxr-np->box.minr+1.0));
			else 
				d=(DenMatrix[nbox.maxl-1][nbox.minr-1]-DenMatrix[nbox.maxl-1][nbox.minr-2]-DenMatrix[nbox.minl-2][nbox.maxr-1]+DenMatrix[nbox.minl-2][nbox.maxr-2])/((np->box.maxl-np->box.minl+1)*(np->box.maxr-np->box.minr+1.0));
			if(d<=0.8&&np->ObCount>=0.8*v*v/pow(4.0,splittimes)&&h<splittimes)
				SplitQuadNode(pt,np,splittimes,AdjMatrix);
			else if(np->ObCount<0.8*v*v/pow(4.0,splittimes))
			{
				free(np);
				p->children[i]=NULL;
			}
		}
	}
	else 
	{
		for(int i=0;i<4;++i)
			p->children[i]=NULL;
	}
	return;
}




void InitQuadNode(ObInf splitpoint,QuadNode *p)
{
	QuadBox nbox;
	nbox.minl=p->box.minl;
	nbox.minr=p->box.minr;
	nbox.maxl=splitpoint.l;
	nbox.maxr=splitpoint.r;
	p->children[2]=AddNode(p->children[2],nbox,p);
	nbox.minl=p->box.minl;
	nbox.minr=splitpoint.r+1;
	nbox.maxl=splitpoint.l;
	nbox.maxr=p->box.maxr;
	p->children[3]=AddNode(p->children[3],nbox,p);
	nbox.minl=splitpoint.l+1;
	nbox.minr=p->box.minr;
	nbox.maxl=p->box.maxl;
	nbox.maxr=splitpoint.r;
	p->children[1]=AddNode(p->children[1],nbox,p);
	nbox.minl=splitpoint.l+1;
	nbox.minr=splitpoint.r+1;
	nbox.maxl=p->box.maxl;
	nbox.maxr=p->box.maxr;
	p->children[0]=AddNode(p->children[0],nbox,p);
	return;
}



int NoisyCount(//double epsilon,			//隐私预算**************************************拉普拉斯噪声**********************
			   double ranNum,				//随机数
			   int ObCount,					//叶子节点中1的个数
			   int i,						//当前树的高度
			   int h						//树的高度
			   )
{


	//double epsilon=0.5;
	double epsilon=1;

	double delta=1;	//全局敏感度
	epsilon=pow(2,i/3.0)*(pow(2,1/3.0)-1)*epsilon/(pow(2,(h+1)/3.0)-1);

	double b = delta / epsilon;
	//先产生0、1之间的随机数，保留4位小数


	double lapNoise;//分别是概率与噪音量



	//得到噪音值*************************************************************
	//噪音量，符合拉普拉斯分布

	if(ranNum < 0.5)//噪音量为负
		lapNoise = log(2 * ranNum) * b;
	else //噪音量为正
		lapNoise = -log(2 * (1 - ranNum)) * b;

	//***********************************************************************
	cout<<"添加的噪音量："<<lapNoise<<endl;
	
	int NoisyCount=ObCount+(int)lapNoise;
	return NoisyCount;
}
double RanNumber(int n)															//产生随机数
{

	double ranNum;//0到1之间均匀分布的随机数，即为拉普拉斯分布值
	srand( (unsigned)time( NULL ) );//srand()函数产生一个以当前时间开始的随机种子 


    //for(int k=0;k<50;k++)
	cout<<"n:"<<n<<endl;
	for(int k=0;k<n;k++)
		ranNum = rand()%1000;
	ranNum=ranNum / 1000.0;//0到1之间均匀分布的随机数

	cout<<"随机数的值："<<ranNum<<endl;
	return ranNum;

}




vector<double> BubbleSort(vector<double> p,int n)//冒泡排序
{
	int i ,j;
	//bool ischanged;//设计跳出条件
	for(j=n;j>0;j--)
	{
		//ischanged =false;
		for(i=0;i<j-1;i++)
		{
			if(p[i]>p[i+1])//如果发现较重元素就向后移动
			{
				double temp=p[i];
				p[i]=p[i+1];
				p[i+1]=temp;
				//ischanged =true;
			}
		}
		//if(!ischanged)//若没有移动则说明序列已经有序，直接跳出
			//break;
	}
	return p;
}


//指数机制
int ExpCount(int i,//当前树的高度
			 int h,//树的高度
			 int v,//矩阵大小
			 int n,//分裂点个数，分裂点数组的大小
			 vector<double> dis//存储分裂点对应密度差
			 )
{
	double epsilon=1;//隐私预算


	//敏感度函数
	double delta;
	delta=pow(4.0,i+1)/pow((double)v,2);

	double sum=0;//将放大后的概率相加，使最后的概率和为1

	vector<double> p;
	//vector<double> temp;
	//int k=0;//dis[]的大小



	for(int k=0;k<n;k++)//算出概率值、概率和
	{
		p.push_back(exp(epsilon*delta/(2*h*delta)*dis[k]));
		sum+=p[k];
	}

	for(int k = 0; k < n; k++)//算出分布值
	{
		p[k] /= sum;
		
	}
	vector<double> temp(p);
	//产生随机数
	double ranNum=RanNumber(n+20);
	
	//将概率值排序
	//quick_sort(p,0,n-1);
	vector<double> p1(BubbleSort(p,n));
	int pm=0;

	for(int k = 1; k < n-1; k++)//算出分布值
	{
		if(ranNum<p1[k-1])
		{
			//cout<<p[k]<<endl;
			//cout<<k+1<<endl;//a[]中对应分裂点坐标
			pm=k;
			break;
		}
		p1[k]+=p1[k-1];
	}
	for(int k = 0; k < n; k++)//算出分布值
	{
		if(temp[k]==p1[pm])
		{
			pm=k;
			break;
		}

	}
	cout<<pm<<endl;
	return pm;
}