#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <set>
#include <algorithm>
#include <sstream>
using namespace std;

map< string, vector<int> > RefIndex;//index for ref
map<string, vector<int> >::iterator ItRefIndex;
string fa_whole;//string of ref.fa
string ref_name;

char CHR1[500], CHR2[500];
int CHR1_begin, CHR1_end, CHR2_begin, CHR2_end;
char ref_ori_1, ref_ori_2;//orientation of reference gene

int totalMis;
int halfLength;
int numSpPair;//number of supporting 
map<string,string> InputFastq; //input reads
map<string, vector<int> > InputFastqPart;//patial information    part(1,2) orien(1,2) position(int)
map<string, vector<int> >::iterator ItFastqPart;
map<string,string> InputQ;//quality

map<string, vector<vector<int> > > Gapped;//orien(1,2) start alignment1 gap alignment2
map<string, vector<vector<int> > >::iterator ItGapped;

vector<int> PosiPair;//list of mapped region
set<int> Annotation1;//list of annotated exon from ucsc
set<int> Annotation2;
int tile_length;

pair<int,int> PositionPair;
map<pair<int,int>, vector<string> > PositionFilter;
map<pair<int,int>, vector<string> >::iterator ItPositionFilter;

set<vector<int> > PCR;

set<vector<int> >::iterator ItPCR;
int MINOVLP;

void usage() {
        cout<<
                "singleGappedAlignment - find gapped alignment and output CIGAR\n"
                "    ";
}

string revString(string Forward){
        string Reverse;
        for(int i=Forward.length()-1;i>=0;i--){
                if(Forward[i]=='A')
                        Reverse.push_back('T');
                else if(Forward[i]=='T')
                        Reverse.push_back('A');
                else if(Forward[i]=='G')
                        Reverse.push_back('C');
                else if(Forward[i]=='C')
                        Reverse.push_back('G');
        }
        return (Reverse);
}

string revQ(string A){
	string B;
	for(int i=A.length()-1;i>=0;i--)
		B.push_back(A[i]);
	return (B);
}

//index the merged fa file with length 6
void getRefIndex(ifstream& input_ref, int mask){
	string fa_line;
	getline(input_ref, fa_line);
	ref_name = fa_line.substr(1);
	ref_ori_1 = fa_line[fa_line.find('[')+1];
	ref_ori_2 = fa_line[fa_line.find('[')+2];
	sscanf(ref_name.c_str(), "%[^:]:%d-%d+%[^:]:%d-%d[%*c%*c](%d)", CHR1, &CHR1_begin,&CHR1_end, CHR2, &CHR2_begin, &CHR2_end,&numSpPair);
	while(!input_ref.eof()){
		getline(input_ref, fa_line);
		fa_whole += fa_line;
	}
	if(mask != 1) transform(fa_whole.begin(),fa_whole.end(),fa_whole.begin(),::toupper);
	for(int i=0;i<int(fa_whole.length());i++){
		RefIndex[ fa_whole.substr(i, tile_length) ].push_back(i);
	}
}
		
int HammingDis(string A, string B, int mis){
	if(A.length()!=B.length()) cout << "Hamming wrong  " <<A<<'\t'<<B<<endl;
	int Dis=0;
	for(int i=0; i<int(A.length()); i++){
		if(A[i]!=B[i])
			Dis++;
		if(Dis > mis)
			return (100);
	}
	return(Dis);
}
			


void extendSegment1(string Name, string singleRead, int start, int segLength, int ori, int mismatch, int segMis){
	int flag = 0;
	int Hamm1, Hamm2;
//	int tempHamm;
	if(ori == 2){
		singleRead=revString(singleRead);
	}
	int Len = singleRead.length();
	int Len_ref = fa_whole.length();
	vector<int> information;
	for(int i=0; i<=Len-segLength-MINOVLP;i++){
		if(start+segLength+i+2 > Len_ref) continue;
		string canonical1 = fa_whole.substr(start+segLength+i, 2);
		if(canonical1 != "GT" && canonical1 != "CT")
			continue;
		string s = singleRead.substr(segLength-tile_length+i,tile_length);
		vector<int>::iterator it = RefIndex[s].begin();
		while(it!=RefIndex[s].end()){
			if(*it == start+segLength-tile_length+i){
				
				if(start+segLength-tile_length+i > Len_ref){
					flag = 1;
					break;
				}
				if(i <= tile_length)
					Hamm1 = segMis;
				else
					Hamm1 = HammingDis(singleRead.substr(0+segLength, i-tile_length), fa_whole.substr(start+segLength, i-tile_length),mismatch) + segMis;
				//	Hamm1 += HammingDis(singleRead.substr(0+segLength+i-tile_length-1, 1), fa_whole.substr(start+segLength+i-tile_length-1, 1), mismatch);
				if(Hamm1 <= mismatch){
					string s2 = singleRead.substr(segLength+i,tile_length);
					vector<int>::iterator it2 =  RefIndex[s2].begin();
					while(it2 != RefIndex[s2].end()){
						if(*it2 <= *it){it2++; continue;}
						if(*it2+tile_length+Len-(segLength+i+tile_length) > Len_ref){it2++; continue;}
						string canonical2 = fa_whole.substr (*it2-2,2);
						if( (canonical1 == "GT"&& canonical2 == "AG") || (canonical1 == "CT"&& canonical2 == "AC")){
						Hamm2 = HammingDis(singleRead.substr(segLength+i+tile_length,Len-(segLength+i+tile_length)), fa_whole.substr(*it2+tile_length, Len-(segLength+i+tile_length)),mismatch-Hamm1);
						if(Hamm2 <= mismatch-Hamm1){
							information.clear();
							information.push_back(ori);
							information.push_back(start);
							information.push_back(segLength+i);//al1
							information.push_back(*it2 - *it -tile_length);//gap
							information.push_back(Len-segLength-i);//al2
							information.push_back(Hamm1);
							information.push_back(Hamm2);
							Gapped[Name].push_back(information);
						}
						}
						it2++;
					}
				}
				else flag = 1;
				break;
			}
			it++;
		}
		if(flag == 1) break;
	}
}
				
void extendSegment2(string Name, string singleRead, int start, int segLength, char ori, int mismatch, int segMis){
	int Hamm1, Hamm2;
        int flag = 0;
	if(ori == 2) singleRead=revString(singleRead);
	int Len = singleRead.length();
        vector<int> information;
        for(int i=0; i<=Len-segLength-MINOVLP;i++){
		if(start-i-2 < 0) continue;
		string canonical1 = fa_whole.substr(start-i-2, 2);
                if(canonical1 != "AG" && canonical1 != "AC")
                        continue;
                string s = singleRead.substr(Len-segLength-i,tile_length);
                vector<int>::iterator it = RefIndex[s].begin();
                while(it!=RefIndex[s].end()){
                        if(*it == start-i){
				
				if(start-i+tile_length<0){
                                        flag = 1;
                                        break;
                                }
				if(i <= tile_length)
                                        Hamm1 = segMis;
                                else
                                       Hamm1 = HammingDis(singleRead.substr(Len-segLength-i+tile_length, i-tile_length), fa_whole.substr(start-i+tile_length, i-tile_length), mismatch) + segMis;
                                if(Hamm1 <= mismatch){
                                        string s2 = singleRead.substr(Len-segLength-i-tile_length,tile_length);
                                        vector<int>::iterator it2 =  RefIndex[s2].begin();
                                        while(it2 != RefIndex[s2].end()){
                                                if(*it2 >= *it){it2++; continue;}
						if(*it2-(Len-segLength-i-tile_length)<0){it2++; continue;}
						string canonical2 = fa_whole.substr (*it2+tile_length,2);
                                                if( (canonical1 == "AG"&& canonical2 == "GT") || (canonical1 == "AC"&& canonical2 == "CT")){
                                                Hamm2 = HammingDis(singleRead.substr(0,Len-segLength-i-tile_length), fa_whole.substr(*it2-(Len-segLength-i-tile_length), Len-segLength-i-tile_length), mismatch-Hamm1); 
                                                if(Hamm2 <= mismatch-Hamm1){
                                                        information.clear();
                                                        information.push_back(ori);
                                                        information.push_back(*it2-(Len-segLength-i-tile_length));
                                                        information.push_back(Len-segLength-i);//al1
                                                        information.push_back(*it - *it2 -tile_length);//gap
                                                        information.push_back(segLength+i);//al2
                                                        information.push_back(Hamm2);
                                                        information.push_back(Hamm1);
                                                        Gapped[Name].push_back(information);
                                                     
                                                }
						}
                                                it2++;
                                        }
                                }
                       		else flag = 1;
				break;
                        }
                        it++;
                }
		if(flag == 1) break;
        }
}
			
void readFq(ifstream& input_fq){
	string fq_line;
	string Name_fq;
	while(!input_fq.eof()){
		getline(input_fq, fq_line);
		if(fq_line.length()==0)  //in case of a blank line
			break;
		Name_fq = fq_line.substr(1);
		getline(input_fq, fq_line);
		InputFastq[Name_fq] = fq_line;
	
		getline(input_fq, fq_line);
		getline(input_fq, fq_line);
		InputQ[Name_fq] = fq_line;
	}
}


void readBwt(ifstream& input_bwt){
	string bwt_line;
	string Name_Index;
	char Orien;
	int Position;
	while(!input_bwt.eof()){
		getline(input_bwt, bwt_line);
		if(bwt_line.length()==0)//in case of blank line
			break;
		istringstream ss(bwt_line);
		string Chr,temp1,temp2,err;
		ss >> Name_Index >> Orien >> Chr >> Position >> temp1 >> temp2 >> err;

		int part = atoi( Name_Index.substr(Name_Index.length()-1).c_str() );
		string Name = Name_Index.substr(0,Name_Index.length()-2);
		InputFastqPart[Name].push_back(part);
		int Ori;
		if(Orien == '+'){Ori=1;}
		if(Orien == '-'){Ori=2;}
		InputFastqPart[Name].push_back(Ori);
		InputFastqPart[Name].push_back(Position);
		int segMis = int(err.length()/5);
		InputFastqPart[Name].push_back(segMis);
	}
}


void GappedAlign(){
	for(ItFastqPart = InputFastqPart.begin();ItFastqPart != InputFastqPart.end();ItFastqPart++){
		string Name = ItFastqPart->first;
		for(int i=0;i*4!=int(ItFastqPart->second.size());i++){
			int part = ItFastqPart->second[i*4];
			int ori = ItFastqPart->second[i*4+1];
			int position = ItFastqPart->second[i*4+2];
			int segMis = ItFastqPart->second[i*4+3];
			if((part==1&&ori==1)||(part==2&&ori==2)){
				extendSegment1(Name,InputFastq[Name],position,halfLength,ori,totalMis,segMis);
			}
			else{
				extendSegment2(Name,InputFastq[Name],position,halfLength,ori,totalMis,segMis);
			}
				
		}
	}
}
	
void Filter_Canonical(){
	vector< vector<int> >::iterator it;
        for(ItGapped=Gapped.begin();ItGapped!=Gapped.end();){
        
		
		
                for(it = ItGapped->second.begin(); it != ItGapped->second.end(); ){
			
                       if((fa_whole.substr( (*it)[1]+(*it)[2],2 ) == "GT" && fa_whole.substr( (*it)[1]+(*it)[2]+(*it)[3]-2,2 ) == "AG") || (fa_whole.substr( (*it)[1]+(*it)[2],2 ) == "CT" && fa_whole.substr( (*it)[1]+(*it)[2]+(*it)[3]-2,2 ) == "AC"))

{
it++;
continue;

}
			
			ItGapped->second.erase(it);
                }
		if(ItGapped->second.size() == 0)
			Gapped.erase(ItGapped++);
		else ItGapped++;
        }

}
void Get_ExonPred(ifstream& input_fullbwt){
	string line;
	vector<int> Posi;
	while(!input_fullbwt.eof()){
		getline(input_fullbwt, line);
		if(line.length()==0)//in case of blank line
                        break;
                istringstream ss(line);
                string a,b,c;
		int Position;
                ss >> a >> b >> c >> Position;
		Posi.push_back(Position);
	}
	sort(Posi.begin(),Posi.end());
	int flag = 0;
	for(int i = 0; i <int(Posi.size()); i++){
		if(flag == 0){
			PosiPair.push_back(Posi[i]);
			if(i==int(Posi.size()-1)){
				PosiPair.push_back(Posi[i]+50);
	
				break;
			}	
			flag=1;
		}
		if(i==int(Posi.size()-1)){
			PosiPair.push_back(Posi[i]+50);
                                break;
		}
		if(Posi[i]>Posi[i+1]-300)
			continue;
		else{
			PosiPair.push_back(Posi[i]+50);
			flag=0;
		}
	}
}

void Get_ExonAnno(ifstream& input_genelist){
	string line;
	while(!input_genelist.eof()){
		getline(input_genelist, line);
		if(line.length()==0)
			break;
		istringstream ss(line);
		string a;
		string chr,chr1,chr2;
		chr1=CHR1;
		chr2=CHR2;
		int begin, end;
		ss >> a;

		if(a == "@"){
			ss >> chr >> begin >> end;
			if(chr==chr1){
				if(begin >= CHR1_begin && end <= CHR1_end){
					if(ref_ori_1 == '+'){
						Annotation1.insert(1+begin-CHR1_begin);
						Annotation1.insert(end-CHR1_begin);
					}
					else{
						Annotation1.insert(CHR1_end-begin-1);
                                                Annotation1.insert(CHR1_end-end);
					}
				}
			}
			if(chr==chr2){
				if(begin >= CHR2_begin && end <= CHR2_end){
					if(ref_ori_2 == '+'){
					Annotation2.insert(1+begin-CHR2_begin+CHR1_end-CHR1_begin+1);
                                        Annotation2.insert(end-CHR2_begin+CHR1_end-CHR1_begin+1);
					}
					else{
						Annotation2.insert(CHR2_end+CHR1_end-CHR1_begin+1-begin-1);
	                                        Annotation2.insert(CHR2_end+CHR1_end-CHR1_begin+1-end);
					}
				}
			}
			
		else
			continue;
		}
		
	}
	//sort(AnnotationPair.begin(),AnnotationPair.end());
	//AnnotationPair.erase(unique(AnnotationPair.begin(),AnnotationPair.end()),AnnotationPair.end());
}
	
/*void Filter_ExonAnnoPred(){
	vector< vector<int> >::iterator it;
        for(ItGapped=Gapped.begin();ItGapped!=Gapped.end();ItGapped++){
		int flag = 0;
	        int flag2 = 0;

                for(it = ItGapped->second.begin(); it != ItGapped->second.end();it++ ){
			for(int i = 0; i< int(AnnotationPair.size());i++){
				if(flag==0){
					if((*it)[1] < AnnotationPair[i]){
						if(i%2 == 0)
							break;
						else
							flag=1;
					}
				}
				if(flag==1){
					if((*it)[1]+(*it)[2]+(*it)[3] < AnnotationPair[i]){

                                                if(i%2 == 0)
                                                        break;
                                                else{
                                                        flag=2;
                                                        break;
                                                }
                                        }
                                }
			}
			for(int i = 0; i< int(PosiPair.size());i++){
                                if(flag2==0){
                                        if((*it)[1] < PosiPair[i]){

                                                if(i%2 == 0)
                                                        break;
                                                else
                                                        flag2=1;
                                        }
                                }
                                if(flag2==1){
                                        if((*it)[1]+(*it)[2]+(*it)[3] < PosiPair[i]){
                                                if(i%2 == 0)
                                                        break;
                                                else{
                                                        flag2=2;
                                                        break;
                                                }
                                        }
                                }
                        }
			(*it).push_back(flag);
			(*it).push_back(flag2);
		}
	}
}

*/		
void Filter_Final(){
	 vector< vector<int> >::iterator it;	
	
	for(ItGapped=Gapped.begin();ItGapped!=Gapped.end();){
                for(it = ItGapped->second.begin(); it != ItGapped->second.end();){
			if((*it)[2]<MINOVLP||(*it)[4]<MINOVLP)
				ItGapped->second.erase(it);
//			else if(((*it)[2]<14&&(*it)[5]>1)||((*it)[4]<14&&(*it)[6]>1))
//				ItGapped->second.erase(it);
			else it++;
		}
		if(ItGapped->second.size() == 0){
                        Gapped.erase(ItGapped++);
			continue;
		}
		if(ItGapped->second.size() > 1){
			int p;
          	        int min = 100;
			int c = 0;
			for(it = ItGapped->second.begin(); it!=ItGapped->second.end();it++,c++ ){
                                int mis = (*it)[5]+(*it)[6];
                                if(mis<=min) { min=mis;p=c;}
        	        }
			c = 0;
			for(it = ItGapped->second.begin(); it!=ItGapped->second.end();c++ ){
                                if(c != p)
					ItGapped->second.erase(it);
				else it++;
               	        }
		}
		vector <int> temp;
		temp.push_back(ItGapped->second[0][1]);
		temp.push_back(ItGapped->second[0][2]);
		temp.push_back(ItGapped->second[0][3]);
		temp.push_back(ItGapped->second[0][4]);
		ItPCR = PCR.find(temp);
		if(ItPCR != PCR.end()){
			Gapped.erase(ItGapped++);
                        continue;
		}
		else PCR.insert(temp);
			
		PositionPair.first = ItGapped->second[0][1]+ItGapped->second[0][2];
                PositionPair.second = ItGapped->second[0][3];
                PositionFilter[PositionPair].push_back(ItGapped->first);
                ItGapped++;
	}
	for(ItPositionFilter=PositionFilter.begin();ItPositionFilter!=PositionFilter.end();ItPositionFilter++){
/*		if(ItPositionFilter->second.size() == 1){
		//	if((Gapped[ItPositionFilter->second[0]][0][2]<15&&Gapped[ItPositionFilter->second[0]][0][5]>0) || (Gapped[ItPositionFilter->second[0]][0][4]<15&&Gapped[ItPositionFilter->second[0]][0][6]>0))
			//if((Gapped[ItPositionFilter->second[0]][0][2]<15) || (Gapped[ItPositionFilter->second[0]][0][4]<15))
			//Gapped.erase(ItPositionFilter->second[0]);
			int p[2];
			p[0] = Gapped[ ItPositionFilter->second[0] ][0][1] + Gapped[ ItPositionFilter->second[0] ][0][2]-1;
			p[1] = Gapped[ ItPositionFilter->second[0] ][0][3] + p[0]+1;
			if( Annotation1.find(p[0]) == Annotation1.end() || Annotation2.find(p[1]) == Annotation2.end() ){
//	if(ItPositionFilter->second[0] == "SRR018261.3940550/1")
//		cout << p[0] <<'\t'<<p[1]<<'\t'<<*Annotation1.find(p[0])<<'\t'<<*Annotation2.find(p[1])<<endl;
				Gapped.erase(ItPositionFilter->second[0]);
			}
		}
*/
		if(ItPositionFilter->second.size() < 8){
			int p[2];
                        p[0] = Gapped[ ItPositionFilter->second[0] ][0][1] + Gapped[ ItPositionFilter->second[0] ][0][2]-1;
                        p[1] = Gapped[ ItPositionFilter->second[0] ][0][3] + p[0]+1;
			if(ItPositionFilter->second.size()>1){
			if( Annotation1.find(p[0]) == Annotation1.end() && Annotation2.find(p[1]) == Annotation2.end() ){
				for(int i=0;i<int(ItPositionFilter->second.size());i++){
					Gapped.erase(ItPositionFilter->second[i]);
				}
				continue;
			}
			if( (Annotation1.find(p[0]) == Annotation1.end() || Annotation2.find(p[1]) == Annotation2.end())){
				if(numSpPair<5 || ItPositionFilter->second.size()<3){
					for(int i=0;i<int(ItPositionFilter->second.size());i++){
                                        Gapped.erase(ItPositionFilter->second[i]);
                                	}
                                	continue;
				}
				int max1=0,max2=0;
				for(int i=0;i<int(ItPositionFilter->second.size());i++){
					max1 = max1 > Gapped[ ItPositionFilter->second[i] ][0][2] ? max1:Gapped[ ItPositionFilter->second[i] ][0][2];
					max2 = max2 > Gapped[ ItPositionFilter->second[i] ][0][4] ? max2:Gapped[ ItPositionFilter->second[i] ][0][4];
                                }
				if(max1<15 || max2 <15 ){
					for(int i=0;i<int(ItPositionFilter->second.size());i++){
        	                                Gapped.erase(ItPositionFilter->second[i]);
                	                }
                        	        continue;
				}
			}
			}
			if(ItPositionFilter->second.size()==1){
			if( Annotation1.find(p[0]) == Annotation1.end() || Annotation2.find(p[1]) == Annotation2.end() ){
                                
                                        Gapped.erase(ItPositionFilter->second[0]);
                                
                        }
                        }
		}
	}

}
int ifAnnotated(int leftPosition, int rightPosition){
	int l,r;
	if(Annotation1.find(leftPosition) == Annotation1.end())
		l = 0;
	else l = 1;
	if(Annotation2.find(rightPosition) == Annotation2.end())
                r = 0;
        else r = 1;
	if(l == 0 && r == 0)return 0;
	else if(l == 0 && r == 1)return 2;
	else if (l == 1 && r == 0)return 1;
	else return 3;
}
						



void out2Sam(){
	vector< vector<int> >::iterator it;
        for(ItGapped=Gapped.begin();ItGapped!=Gapped.end();ItGapped++){
                cout<<ItGapped->first;
		if(ItGapped->second.size()==1){
			if(ItGapped->second[0][0]==1)
				cout << "\t0\t"<<ref_name<<'\t'<<ItGapped->second[0][1]+1<<"\t255\t"<< ItGapped->second[0][2]<<"M"<<ItGapped->second[0][3]<<"N"<<ItGapped->second[0][4]<<"M"<<"\t*\t0\t0\t"
					<<InputFastq[ItGapped->first]<<'\t'<<InputQ[ItGapped->first]<<'\t'<<"NM:i:0\tXX:A:"
					<<ifAnnotated( ItGapped->second[0][1]+ItGapped->second[0][2]-1, ItGapped->second[0][1]+ItGapped->second[0][2]+ItGapped->second[0][3] )<<"\tNH:i:1"<<endl;
			if(ItGapped->second[0][0]==2)
                                cout << "\t16\t"<<ref_name<<'\t'<<ItGapped->second[0][1]+1<<"\t255\t"<< ItGapped->second[0][2]<<"M"<<ItGapped->second[0][3]<<"N"<<ItGapped->second[0][4]<<"M"<<"\t*\t0\t0\t"
                                        <<revString(InputFastq[ItGapped->first])<<'\t'<<revQ(InputQ[ItGapped->first])<<'\t'<<"NM:i:0\tXX:A:"
					<<ifAnnotated( ItGapped->second[0][1]+ItGapped->second[0][2]-1, ItGapped->second[0][1]+ItGapped->second[0][2]+ItGapped->second[0][3] )<<"\tNH:i:1"<<endl;
		}
		else{
			int c = 0;int p;
			int min = 100;
               		for(it = ItGapped->second.begin(); it!=ItGapped->second.end();it++,c++ ){
				int mis = (*it)[5]+(*it)[6];
				if(mis<=min) { min=mis;p=c;}
        	        }
			if(ItGapped->second[p][0]==1)
                                cout << "\t0\t"<<ref_name<<'\t'<<ItGapped->second[p][1]+1<<"\t255\t"<< ItGapped->second[p][2]<<"M"<<ItGapped->second[p][3]<<"N"<<ItGapped->second[p][4]<<"M"<<"\t*\t0\t0\t"
                                        <<InputFastq[ItGapped->first]<<'\t'<<InputQ[ItGapped->first]<<'\t'<<"NM:i:0\tXX:A:"
					<<ifAnnotated( ItGapped->second[p][1]+ItGapped->second[p][2]-1, ItGapped->second[p][1]+ItGapped->second[p][2]+ItGapped->second[p][3] )<<"\tNH:i:1"<<endl;
                        if(ItGapped->second[p][0]==2)
                                cout << "\t16\t"<<ref_name<<'\t'<<ItGapped->second[p][1]+1<<"\t255\t"<< ItGapped->second[p][2]<<"M"<<ItGapped->second[p][3]<<"N"<<ItGapped->second[p][4]<<"M"<<"\t*\t0\t0\t"
                                        <<revString(InputFastq[ItGapped->first])<<'\t'<<revQ(InputQ[ItGapped->first])<<'\t'<<"NM:i:0\tXX:A:"
					<<ifAnnotated( ItGapped->second[p][1]+ItGapped->second[p][2]-1, ItGapped->second[p][1]+ItGapped->second[p][2]+ItGapped->second[p][3] )<<"\tNH:i:1"<<endl;
       	 	}
	}
}

int main(int argc, char *argv[]){
	tile_length = atoi(argv[argc-7]);
	MINOVLP = atoi(argv[argc-8]);
	totalMis = atoi(argv[argc-9]);
	halfLength = atoi(argv[argc-10]);
	int mask = atoi(argv[argc-6]);
	ifstream input_ref(argv[argc-2]);
	getRefIndex(input_ref,mask);
	
	ifstream input_fq(argv[argc-3]);
	ifstream input_bwt(argv[argc-1]);
	ifstream input_fullbwt(argv[argc-4]);
	ifstream input_genelist(argv[argc-5]);
	readFq(input_fq);
	readBwt(input_bwt);
	GappedAlign();
	//Get_ExonPred(input_fullbwt);

	Get_ExonAnno(input_genelist);
	//Filter_Canonical();	
//	Filter_ExonAnnoPred();
	Filter_Final();
	out2Sam();
	return (0);

}
	

