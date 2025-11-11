#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>

using namespace std;

class NODE
{
	public:
		int value;
		int group;
		NODE* left, *right;
		NODE(){group=-1; left=right=NULL;}
		NODE(int p) {group=-1; value=p; left=right=NULL;}
};

class GeneInfo
{
	public:
		string geneName;
		int nodeIndex;
		GeneInfo(){geneName=""; nodeIndex=-1; }
};

class TreeAnalysis
{
	public:

		vector<string> Label;
		string speciesTree;
		vector<vector<int> > groups;
		vector<vector<string> > geneName;
		int N, S;

		vector<string> GeneBirth;
		vector<string> GeneDuplication;
		vector<int> GeneLoss;

		int groupIndex;

		void IdentifyOrthoGroup(NODE* root, int curIndex)
		{
			if(curIndex==-1)
			{
				if(root->left==NULL or root->right==NULL) return;
				if(root->value==0)
				{
					IdentifyOrthoGroup(root->left, curIndex);
					IdentifyOrthoGroup(root->right, curIndex);
				}
				else
				{
					groupIndex++;
					groups.push_back(vector<int>());

					IdentifyOrthoGroup(root->left, groupIndex);
					IdentifyOrthoGroup(root->right, groupIndex);
				}
				return;
			}
			else
			{
				if(root->left==NULL or root->right==NULL)
				{
					groups[curIndex].push_back(root->value);
				}
				else
				{
					IdentifyOrthoGroup(root->left, curIndex);
					IdentifyOrthoGroup(root->right, curIndex);
				}
			}
		}

		void findOrthoGroup(string treeLabel)
		{
			int index=0;
			vector<NODE* > stack;
			for(int i=0; i<treeLabel.size(); i++)
			{
				NODE* newNODE=new NODE();

				if(speciesTree[i]!='N')
				{
					newNODE->value=index++;
					stack.push_back(newNODE);
				}
				else
				{
					NODE* right=stack.back(); stack.pop_back();
					NODE* left=stack.back(); stack.pop_back();

					newNODE->left=left;
					newNODE->right=right;
					newNODE->value=treeLabel[i]-'0';

					stack.push_back(newNODE);
				}
			}
			NODE* root=stack[0];
			
			groups.clear();

			groupIndex=-1;

			IdentifyOrthoGroup(root, -1);

		}

		// First Constructor
		TreeAnalysis(string spt, vector<string> labeling)
		{
			speciesTree=spt;
			Label=labeling;

			N=Label.size();
			S=speciesTree.size();
			
			for(int i=0; i<N; i++)
			{
				cout<<"For tree "<<i<<": "<<endl;
				findOrthoGroup(Label[i]);

				for(int j=0; j<groups.size(); j++)
				{
					cout<<"\tGroup "<<j<<": ";
					for(int k=0; k<groups[j].size(); k++)
						cout<<groups[j][k]<<"\t";
					cout<<endl;
				}
			}
		}

		// Second Constructor
		TreeAnalysis(string spt, vector<string> labeling, vector<vector<string> > treeGeneName)
		{
			speciesTree=spt;
			Label=labeling;
			geneName=treeGeneName;

			N=Label.size();
			S=speciesTree.size();
		}


		void printAnalysis()
		{
			for(int i=0; i<N; i++)
			{
				cout<<"For tree "<<i<<": "<<endl;
				findOrthoGroup(Label[i]);

				for(int j=0; j<groups.size(); j++)
				{
					cout<<"\tGroup "<<j<<": ";
					for(int k=0; k<groups[j].size(); k++)
						cout<<geneName[i][groups[j][k]]<<"\t";
					cout<<endl;
				}
			}
		}

		void printOrthoGroups(ofstream& outfile)
		{
			set<set<string> > sst;

			for(int i=0; i<N; i++)
			{
				findOrthoGroup(Label[i]);

				for(int j=0; j<groups.size(); j++)
				{
					set<string> stmp;
					for(int k=0; k<groups[j].size(); k++)
						if(geneName[i][groups[j][k]]!="")
							stmp.insert(geneName[i][groups[j][k]]);
					sst.insert(stmp);
				}
			}

			for(set<set<string> >::iterator it=sst.begin(); it!=sst.end(); it++)
			{
				if((*it).size()<2) continue;
				for(set<string>::iterator j=(*it).begin(); j!=(*it).end(); j++)
					outfile<<*j<<"\t";
				outfile<<endl;
			}

		}

		// Thread-safe version: writes to stringstream buffer instead of file
		void printOrthoGroups_Buffer(stringstream& buffer)
		{
			set<set<string> > sst;

			for(int i=0; i<N; i++)
			{
				findOrthoGroup(Label[i]);

				for(int j=0; j<groups.size(); j++)
				{
					set<string> stmp;
					for(int k=0; k<groups[j].size(); k++)
						if(geneName[i][groups[j][k]]!="")
							stmp.insert(geneName[i][groups[j][k]]);
					sst.insert(stmp);
				}
			}

			for(set<set<string> >::iterator it=sst.begin(); it!=sst.end(); it++)
			{
				if((*it).size()<2) continue;
				for(set<string>::iterator j=(*it).begin(); j!=(*it).end(); j++)
					buffer<<*j<<"\t";
				buffer<<endl;
			}

		}

		void printDetailedAnalysis()
		{
			for(int i=0; i<N; i++)
			{
				cout<<"For tree "<<i<<": "<<endl;
				findOrthoGroup(Label[i]);

				for(int j=0; j<groups.size(); j++)
				{
					cout<<"\tGroup "<<j<<": ";
					for(int k=0; k<groups[j].size(); k++)
						cout<<geneName[i][groups[j][k]]<<"\t";
					cout<<endl;
				}
			}

			// Detailed Analysis
			string unionTree=string(S, '0');
			for(int i=0; i<S; i++)
			{
				for(int j=0; j<N; j++) if(Label[j][i]=='1')
				{
					unionTree[i]='1';
					break;
				}
			}

			int totalGeneDuplication=0;
			int totalGeneLoss=0;
			int totalGeneBirth=0;

			for(int i=0; i<N; i++)
			{
				vector<pair<int,int> > stack;
				int index=0;
				for(int j=0; j<S; j++)
				{
					if(speciesTree[j]!='N')
					{
						stack.push_back(make_pair(index++, j));
					}
					else
					{
						pair<int,int> curpair=stack.back(); stack.pop_back();
						int right1=curpair.first, right2=curpair.second;
						curpair=stack.back(); stack.pop_back();
						int left1=curpair.first, left2=curpair.second;

						if(Label[i][left2]!=Label[i][j])
						{
							if(Label[i][left2]=='0') totalGeneLoss++;
							else if(unionTree[j]=='1') totalGeneDuplication++;
							else totalGeneBirth++;
						}


						if(Label[i][right2]!=Label[i][j])
						{
							if(Label[i][right2]=='0') totalGeneLoss++;
							else if(unionTree[j]=='1') totalGeneDuplication++;
							else totalGeneBirth++;
						}

						if(left1>=0)
						{
							if(Label[i][left2]!=Label[i][j])
							{
								if(Label[i][left2]=='0')
								{
									cout<<"Gene loss in species "<<left1<<endl;
								}
								else
								{
									if(unionTree[j]=='1')
									{
										cout<<geneName[i][left1]<<" is a duplicated gene."<<endl;
									}
									else
									{
										cout<<geneName[i][left1]<<" is a new created gene."<<endl;
									}
								}
							}
						}
	

						if(right1>=0)
						{
							if(Label[i][right2]!=Label[i][j])
							{
								if(Label[i][right2]=='0')
								{
									cout<<"Gene loss in species "<<right1<<endl;
								}
								else
								{
									if(unionTree[j]=='1')
									{
										cout<<geneName[i][right1]<<" is a duplicated gene."<<endl;
									}
									else
									{
										cout<<geneName[i][right1]<<" is a new created gene."<<endl;
									}
								}
							}
						}

						stack.push_back(make_pair(-1, j));
					}
				}

			}
			cout<<"Total Gene Birth: "<<totalGeneBirth<<endl;
			cout<<"Total Gene Loss: "<<totalGeneLoss<<endl;
			cout<<"Total Gene Duplication: "<<totalGeneDuplication<<endl;
	
		}

		void printGeneInfo()
		{
			// Gene birth, duplication and loss analysis
			string unionTree=string(S, '0');
			for(int i=0; i<S; i++)
			{
				for(int j=0; j<N; j++) if(Label[j][i]=='1')
				{
					unionTree[i]='1';
					break;
				}
			}

			for(int i=0; i<N; i++)
			{
				vector<GeneInfo* > stack;
				int index=0;
				for(int j=0; j<S; j++)
				{
					GeneInfo* newNode=new GeneInfo();
					newNode->nodeIndex=j;

					if(speciesTree[j]!='N')
					{
						newNode->geneName=geneName[i][index++];
					}
					else
					{
						GeneInfo *right=stack.back(); stack.pop_back();
						GeneInfo *left=stack.back(); stack.pop_back();

						if(Label[i][left->nodeIndex]=='1' and Label[i][right->nodeIndex]=='1')
						{
							/*
							stringstream ss;
							string tmp=left->geneName;
							int h=1;
							while( !isalpha(tmp[h])) h++;
							ss<<"S"<<(int)(speciesTree2[j]-'a')<<tmp.substr(h);
							newNode->geneName=ss.str();
							*/
						}
						else if(Label[i][left->nodeIndex]=='0' and Label[i][right->nodeIndex]=='1')
						{
							if(Label[i][j]=='0')
							{
								if(unionTree[j]=='0')
								{
									if(right->geneName!="") GeneBirth.push_back(right->geneName);
									unionTree[j]='1';
								}
								else
								{
									if(right->geneName!="")	GeneDuplication.push_back(right->geneName);
								}
							}
							else
							{
								/*
								string tmp=right->geneName;
								int h=1;
								while( !isalpha(tmp[h])) h++;
								stringstream ss, ss2;
								ss<<"S"<<(int)(speciesTree2[j]-'a')<<tmp.substr(h);
								newNode->geneName=ss.str();
								ss2<<"S"<<(int)(speciesTree2[left->nodeIndex]-'a')<<tmp.substr(h);
								left->geneName=ss2.str();
								*/

								GeneLoss.push_back(left->nodeIndex);
							}
						}
						else if(Label[i][left->nodeIndex]=='1' and Label[i][right->nodeIndex]=='0')
						{
							if(Label[i][j]=='0')
							{
								if(unionTree[j]=='0')
								{
									if(left->geneName!="") GeneBirth.push_back(left->geneName);
									unionTree[j]='1';
								}
								else
								{
									if(left->geneName!="") GeneDuplication.push_back(left->geneName);
								}
							}
							else
							{
								/*
								string tmp=left->geneName;
								int h=1;
								while(!isalpha(tmp[h])) h++;
								stringstream ss, ss2;
								ss<<"S"<<(int)(speciesTree2[j]-'a')<<tmp.substr(h);
								newNode->geneName=ss.str();
								ss2<<"S"<<(int)(speciesTree2[right->nodeIndex]-'a')<<tmp.substr(h);
								right->geneName=ss2.str();
								*/

								GeneLoss.push_back(right->nodeIndex);
							}
						}
					}
					stack.push_back(newNode);
				}
			}
		}
};
