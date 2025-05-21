#ifndef _HUNGARIAN_
#define _HUNGARIAN_

#include <list>
#include <string>
#include <vector>
#include <queue>
#include <iostream>
#include <limits>


using namespace std;


class Hungarian
{
public:
  int n;
  vector<int> matchingX;
  vector<int> matchingY;
  vector<int>Lx;
  vector<int>Ly;
  vector<int>S;
  vector<int>T;
  vector<int>Nl_S;
  list<int>augmentpath;
  int totalweight;
  typedef vector<int> row;
  typedef vector<row> matrix;
  matrix weight;
  matrix El;

  Hungarian ( ) 
  {
    totalweight=0;
    n = 0;
  }

  Hungarian( matrix w )
  { 
    totalweight = 0;
    n = (int)w.size( );
    Lx.clear();
    Ly.clear();
    El.clear();
    S.clear();
    T.clear();
    Nl_S.clear();
 
    //	cout<<"Get weight matrix from parameter"<<endl;
    for (int i=0; i<n; i++)
      weight.push_back(w[i]);
    

    //	cout<<"Generate initial labelling L"<<endl;
    for (int i=0; i<n; i++)
      El.push_back(vector<int>(n,0));
    for (int i=0; i<n; i++)
      Ly.push_back( 0 );
    for (int i=0; i<n; i++)
    {
      int maxweight = 0;
      for (int j=0; j<n; j++)
      {
	if (weight[i][j] > maxweight)
	  maxweight = weight[i][j];
      }	 
      Lx.push_back(maxweight);
    }
    
    //	cout<<"Generate El, where El[x][y]=1 if Lx[x]+Ly[y] = weight[x][y], otherwise 0"<<endl;
    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
      {
	if (Lx[i]+Ly[j] == weight[i][j])
	  El[i][j] = 1;
      }
    }

    //	cout<<"Find the initial maximal matching in El"<<endl;
    for (int i=0; i<n; i++)
    {
      matchingX.push_back(-1);
      matchingY.push_back(-1);
    }
    
    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
      {
	if (El[i][j] == 1 && matchingX[i] == -1 && matchingY[j] == -1)
	{
	  matchingX[i] = j;
	  matchingY[j] = i;
	  break;
	}
      }
    }

    //	cout<<"Main loop"<<endl;

    while(1)
    {
      int v = findFreevertex();
      if (v == -1)
	break; 
      
      //	cout<<"\t If mathcing is not perfect, pick a free vertex u from X, set S<-{u}, T<-empty"<<endl;
      S.clear();
      T.clear();
      Nl_S.clear();
      S.push_back (v);
      for (int i=0; i<n; i++)
      {
	if (El[v][i] ==1)
          Nl_S.push_back (i);
      }

      //	cout<<"\t Nl_S==T"<<endl;
      if (isequal(Nl_S, T))
	updateLabels( );

      //	cout<<"\t Nl_S !=T"<<endl;
      while (!isequal(Nl_S, T))
      {
	int y = diff(Nl_S, T);
        if (isFreeinY(y))
        {
	  augment(v, y);
	  break;
        }
        else
        {
	  int z = matchingY[y];
	  uniqueInsert(S,z);
          uniqueInsert(T,y);
          for (int i=0; i<n; i++)
          {
            if (El[z][i] ==1)
	      uniqueInsert(Nl_S,i);
          }       
	  if (isequal(Nl_S,T))
            updateLabels( );
	}
      }
    }

    //	cout<<"Matching is perfect, calculate total weight"<<endl;
    for (int i=0; i<n; i++)
      totalweight += weight[i][matchingX[i]];    
  }


  bool isFreeinY (int i)
  {
     if(matchingY[i] == -1 )
        return true;
     else
        return false;

  }

  bool isFreeinX (int i)
  {
     if(matchingX[i] == -1 )
        return true;
     else
        return false;
  }

  int findFreevertex( )
  {
    for (int i=0; i<n; i++)
    {
      if (isFreeinX (i))
	return i;
    }
    return -1;
  }
 
  void updateLabels()
  {
    //	cout<<"\t\t Update Labels"<<endl;
    int sizeS = (int)S.size();
    int sizeT = (int)T.size();

    //calculate nT = ~T
    vector<int>nT;
    for (int i=0; i<n; i++)
    {
      int j;
      for ( j=0; j<sizeT; j++)
      {
	if (T[j]==i) 
	  break;
      }
      if (j==sizeT)
	nT.push_back( i );
    }        
    int sizenT = (int)nT.size();

    //calculate alpha
    int alpha = std::numeric_limits<int>::max();
    for (int i=0; i<sizeS; i++)
    {
      int x = S[i];
      for (int j=0; j<sizenT; j++)
      {
	int y = nT[j];
        int tmp = Lx[x] + Ly[y] -weight[x][y];
	if (tmp < alpha)
	  alpha = tmp;	
      }	
    }

    //update labels
    for (int i=0; i<sizeS; i++)
    {
      int x = S[i];
      Lx[x] -= alpha;
    }
    for (int i=0; i<sizeT; i++)
    {
      int y = T[i];
      Ly[y] += alpha;
    }     

    //update El
    for (int i=0; i<n; i++)
    {
      for (int j=0; j<n; j++)
      {
        if(Lx[i]+Ly[j]==weight[i][j])
	  El[i][j] = 1;
	else
	  El[i][j] = 0;
      }
    } 

    //update Nl_S
    Nl_S.clear( );
    for (int i=0; i<sizeS; i++)

    {
      int x = S[i];
      for (int j=0; j<n; j++)
      {
	if (El[x][j] == 1)
	  uniqueInsert(Nl_S,j);
      }
    }
  }

  void augment(int start, int end)
  {
    
    queue<int>Q; //for BFS use
    vector<int>parent(2*n,-1);
 
    //	cout<<"\t\t Find the augment path"<<endl;
    DFS(parent, start, end);
  
    //	cout<<"\t\t Augmenting the matching"<<endl;
    int b = end +n;
    int a = parent[b];
    matchingX[a] = b-n;
    matchingY[b-n] = a;
    while (a!=start)    
    {
      b = parent[a];
      a = parent[b];
      matchingX[a] = b-n;
      matchingY[b-n] = a;       
    }   
  }
  
  void DFS(vector<int>& p, int s, int e)
  {
    vector<int>color(2*n,0);
    DFS_visit( s,e, p, color );
  }

  void DFS_visit(int s, int e,vector<int>& p, vector<int>&color)
  {
    color[s] = 1;
    if(s == e+n)
	return;
    if (s<n)
    {
      for(int i =0; i<n; i++)
      {
	if(El[s][i] == 1)
	{
	  if(color[i+n] == 0)
	  {
	    p[i+n] = s;
	    DFS_visit(i+n,e,p,color);
	  }
	  color[s] = 2;
	}
      }	
    }
    else
    {
      for(int i =0; i<n; i++)
      {
        if(matchingY[s-n] == i)
        {
          if(color[i] == 0)
          {
            p[i] = s;
            DFS_visit(i,e,p,color);
          }
          color[s] = 2;
        }
      }

    }
  }

  int diff(vector<int>a, vector<int>b)
  {
    for (int i=0; i<(int)a.size(); i++)
    {
      int j;
      for (j=0; j<(int)b.size(); j++)
      {
	if (b[j] == a[i])
	  break;
      }
      if (j==(int)b.size())
	return a[i];
    }
    return -1;
  }
   
  bool isequal( vector<int>a, vector<int>b) 
  {
    if (a.size() != b.size())
      return false;
    for (int i =0; i<(int)a.size(); i++)
    {
      int j;
      for(j=0; j<(int)b.size(); j++)
      {
	if( a[i] == b[j] )
	  break;
      }
      if (j==(int)b.size())
	return false;
    }
    return true;
  }

  void uniqueInsert (vector<int> & a, int i)
  {
    for (int j=0;j<(int)a.size();j++)
    {
      if (a[j]==i)
	return;
    }
    a.push_back(i);
  }

  void print()
  {
    for (int i=0; i<(int)matchingX.size(); i++)
    {
      cout<<matchingX[i]<<endl;
    }
    cout<<"Total weight of the optimal mathcing is: "<< totalweight<<endl;
  }
};//endof class Graph  
  
#endif

//Local Variables:
//mode: c++Com
//End:
