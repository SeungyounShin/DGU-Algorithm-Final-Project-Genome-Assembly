#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <map>
#include <algorithm>
#include <tuple>
using namespace std;
//2018112005 컴퓨터공학과 신승윤

string import(){
  string line;
  ifstream file("ref_10000.txt");
  getline (file,line);
  //cout << line << endl;
  file.close();
  return line;
}

void ShortRead(int L, int M,int N) { //shortread 생성
    ofstream fout;
    ifstream fin;
    string shortread, s;
    int idx;
    srand(time(NULL));
    fin.open("ref_10000.txt"); //해당 파일 열거나 없으면 생성
    getline(fin, s); //str 받아오기
    fout.open("shortread_" + to_string(L) + "_" + to_string(M) + ".txt"); //해당 파일 열거나 없으면 생성
    for (int i = 0; i < M; i++) {
        idx = rand() % (s.length() - L); //랜덤 위치 얻기
        shortread = s.substr(idx, L);
        fout << shortread << endl;
        shortread = "";
    }
    fin.close();
    fout.close(); //파일닫기
}

vector<string> importReads(int L, int M, int N){
  vector<string> ret;
  ifstream file("shortread_" + to_string(L) + "_" + to_string(M) + ".txt");

  int count = 0;
  if(file.is_open()){
    while(!file.eof()){
      string tmp;
      if(ret.size() == M) return ret;
      getline(file,tmp);
      ret.push_back(tmp);
      count++;
    }
  }
  file.close();    //파일 닫기
  return ret;
}

int overlap(string s1, string s2,int min_length){
  int start = 0;
  cout << " overlap func :: "<< s1 << " " << s2 << endl;
  start = s1.find(s2.substr(0,min_length));
  cout << start << " "<< s1 << " "<<s2<<endl;
  if(start == -1){
    if(s1.length() >= s2.length()){
      for(int i=1;i<s2.length();i++)
        if(s1.substr(s1.length()-s2.length()+i, s2.length()-i)  == s2.substr(0,s2.length()-i)){
          return s1.length()-s2.length()+i;
        }
    }
    else{
      for(int i=1;i<s1.length();i++)
        if(s1.substr(i,s1.length()-i)  == s2.substr(0,s1.length()-i)){
          return i;
        }
      //is s1 is subset of s2?
    }
    return start;
  }
  return start;
}

int overlap2(string s1, string s2,int min_length){
  int start = -1;
  //case 1 :: full overlap
  start = s1.find(s2);
  if(start > 0 && start < s1.length() && s1.length() > s2.length()){
    return start;
  }
  //case 2 :: ------
  //              --------
  for(int i=1; i<s1.length(); i++){
      string right = s1.substr(i,s1.length()-i);
      string left = s2.substr(0,s1.length()-i);
      //cout << "CASE 2 " << endl;
      //cout << s1 << endl;
      //for(int j=0; j<i;j++) cout << " "; cout << s2 << endl;
      if(right == left){
        return i;
      }
  }

  return start;
}

tuple <int,int,int> pick_maximal_overlap(vector <string> reads, int k, bool prnt){
  bool print = prnt;
  vector <int> out;
  int best_olen = 0;
  int reada, readb;

  //if(print) cout << " pick maximal overlap :: reads " << endl;
  //if(print) for(auto s : reads)cout << " " << s << endl;
  //if(print) cout << "=====" << endl;

  map <string, vector<string> > index;
  int cnt = 0;
  for(auto read : reads){
    for(int i=0; i<read.length()-k+1; i++){
      index[read.substr(i,k)].push_back(read);
      //if(print) cout << read.substr(i,k) << " ";
    }
    cnt ++;
    //if(print) cout << endl;
  }
  //cout<<"flag"<<endl;

  int i,j;i=0;j=0;
  int olenInd = 0;
  string tmpa; string tmpb;
  for(auto r : reads){
    for(auto o : index[r.substr(r.length()-k,k)] ){
      if(print)cout << " r,o :: "<< r << " " << o << endl;
      //if(print)cout << " sub min : " << r.substr(r.length()-k,k) << endl;
      //if(print)for(auto sub : index[r.substr(r.length()-k,k)]) cout << "    " << sub<<endl;
      if(r != o){
        bool flag = false;
        int olen = overlap2(r,o,k);
        int overlap = r.length() - olen;
        if(olen == -1) overlap = 0;
        if(print) cout << " + overlap " << overlap<< " > " << best_olen << " "<< olen << " "<< r << " "<<o<< " " << endl;

        if(overlap > best_olen){
          tmpa = r; tmpb = o;

          if(print) cout << "best overlap " << overlap<< " "<< r << " "<<o<< " " << endl;
          best_olen = overlap;
          olenInd = olen;
        }
      }
    }
  }
  vector <string>::iterator itera = find(reads.begin(), reads.end(), tmpa);
  vector <string>::iterator iterb = find(reads.begin(), reads.end(), tmpb);
  reada = distance(reads.begin(), itera);
  readb = distance(reads.begin(), iterb);
  return make_tuple(reada, readb, olenInd);
}

vector <string> greedySCS(vector <string> shortReads, int k){
  bool print =true;
  clock_t start = clock();
  tuple<int,int,int> out = pick_maximal_overlap(shortReads, k, print);
  int reada = get<0>(out); int readb = get<1>(out); int olen =get<2>(out);
  int count = 0;
  cout << shortReads[reada] << " " << shortReads[readb] << " "<<olen << endl;
  while(olen>-1 && shortReads.size() > 1){
    cout << "progress :: "<<count << " "<< shortReads.size() << endl;
    if(print) cout << " [shortReads-before-erase] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    if(print) cout << "greedySCS - choose "<<  reada << " "<< readb << " "<< olen << endl;
    if(reada == readb) break;
    string read_a = shortReads[reada];
    string read_b = shortReads[readb];
    if(print) cout << "     "<< read_a << " " << read_b << endl;
    if(reada > readb){
      shortReads.erase(shortReads.begin()+reada);
      shortReads.erase(shortReads.begin()+readb);
    }else{
      shortReads.erase(shortReads.begin()+readb);
      shortReads.erase(shortReads.begin()+reada);
    }
    if(print) cout << " [shortReads-after-erase] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    if(read_a.find(read_b) > 0 &&read_a.find(read_b) <= read_a.length() && (read_a.length() > read_b.length())){
      if(print) cout << " push[1] -> " << read_a << " "<<  read_a.find(read_b)<<endl;
      if(print) cout << "       -> " << read_a << " "<<  read_b<< endl;
      shortReads.push_back(read_a);
    }else if(read_b.find(read_a) > 0 &&read_b.find(read_a) <= read_b.length() && (read_b.length() > read_a.length())){
      if(print) cout << " push[2] -> " << read_b << " "<<  read_b.find(read_a)<< endl;
      if(print) cout << "       -> " << read_a << " "<<  read_b<< endl;
      shortReads.push_back(read_b);
    }else{
      if(print) cout << " push[3] -> " << read_a.substr(0, olen) << " + " << read_b << endl;
      shortReads.push_back(read_a.substr(0, olen) + read_b);
    }
    if(print) cout << " [shortReads-after-insert] " << shortReads.size() <<endl;
    if(print) for(auto t : shortReads) cout << " " << t << endl;
    tuple<int,int,int> tmp = pick_maximal_overlap(shortReads, k, false);
    reada = get<0>(tmp); readb = get<1>(tmp); olen =get<2>(tmp);
    count += 1;
  }
  string ret = "";
  clock_t end = clock();
  cout << "Time spend :: "<< end-start<<endl;
  return shortReads;
}

int main(){
  int M, L,N;
  L = 30;
  M = 400;
  N = 1000;
  string ref = import();
  ShortRead(L, M, N);
  vector<string> shortReads = importReads(L, M, N);

  cout << " # of ShortReads :: "<< shortReads.size() << endl;
  //string ref = "AGATAAATGGGCCCTGCG";
  //string arr[] = {"ACGGATGAGC", "GAGCGGA", "GAGCGAG"};
  //string arr[] = {"AGATAACTG", "TAACTGGGC", "TAACTGGGCCC","AACTGGGC","TGGGCCCTACG"}; // Initialize vector with a string array
  //vector<string> shortReads(arr, arr + sizeof(arr)/sizeof(string));

  vector<string> recon = greedySCS(shortReads, 1);
  cout << " + recon size :: "<<recon.size() << endl;

  ofstream file("ref.txt");
  ofstream file2("recon_vec.txt");
  file << ref;
  for(auto r : recon){
      file2 << r << endl;
  }
  file.close();
  file2.close();

  return 0;
}
