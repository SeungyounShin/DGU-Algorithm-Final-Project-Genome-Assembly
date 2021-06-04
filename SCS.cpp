#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <time.h>
#include <map>
#include <algorithm>
using namespace std;
//2018112005 컴퓨터공학과 신승윤

void init(int length){
  //create RNA seq
  ofstream file;
  file.open("RNA_seq.txt");
  srand(time(NULL));
  for(int i = 0; i<length; i++){
      int randNum = rand() % 4;
      if(randNum == 0) file << 'A';
      else if(randNum == 1) file << 'C';
      else if(randNum == 2) file << 'G';
      else file << 'T';
  }
  file.close();
}

string import(){
  string line;
  ifstream file("RNA_revise.txt");
  getline (file,line);
  file.close();
  return line;
}

void reviseSeq(){
  string line;
  srand(time(NULL));
  ifstream file("RNA_seq.txt");
  getline (file,line);
  for(int i=0; i<line.length()-1; i++){
    if(line[i] == line[i+1]){
      int randNum = rand() % 4;
      char c;
      if(randNum == 0) c='A';
      else if(randNum == 1) c ='C';
      else if(randNum == 2) c='G';
      else c='T';
      line[i] = c;
    }
  }
  file.close();
  ofstream file2;
  file2.open("RNA_revise.txt");
  file2 << line;
  file2.close();
  return ;
}

vector<string> Sequencer(int n, int k, string s){
  vector <string> rtn;
  vector <int> rands;

  random_device rd1,rd2,rd3,rd4;
  mt19937 mt(rd1());
  mt19937 mt2(rd2());
  mt19937 mt3(rd3());
  mt19937 mt4(rd4());
  uniform_int_distribution<int> dist(0, s.length()-k-2);
  uniform_int_distribution<int> dist2(0, 2);
  uniform_int_distribution<int> dist3(0, k-1);
  uniform_int_distribution<int> dist4(0, 4);

  //make short read
  bool check[s.length()];
  char myseq[s.length()];
  for(int i=0; i<s.length(); i++) check[i] = false;
  for(int i=0; i<s.length(); i++) myseq[i] = 'x';
  for(int i=0; i<n; i++){
    int randNum = dist(mt);
    int cover = 0;
    for(int j=0; j<k ;j++) if(check[randNum + j]) cover+=1;
    if(cover >= k-10){i--; continue;}//for uniform distributed reads
    rands.push_back(randNum);
    string sub = s.substr(randNum, k);
    for(int j=0; j<k; j++){
      myseq[randNum+j] = sub[j];
      check[randNum+j] = true;
    }
    //for(int j=0; j<s.length() ;j++) cout << check[j] ;
    //cout << endl;
    rtn.push_back(sub);
  }

  cout << "[shortReads maked] :: "<< rtn.size() <<endl;
  int matched = 0;
  int is_uniform = 0;
  for(int i=0; i<s.length(); i++){
    if(s[i] == myseq[i]) matched += 1;
    if(check[i]) is_uniform+=1;
  }
  cout << " matched :: " << matched << "/" << s.length() << endl;
  cout << " cover   :: " << is_uniform << "/" << s.length() << endl;
  // save to txt
  ofstream file;
  file.open("shortReads.txt");
  // noise (SNIP) to short reads
  for(auto read : rtn){
    int noiseNum = dist2(mt2);
    for(int j=0; j<noiseNum; j++){
      int randInd = dist3(mt3);
      int randChar = dist4(mt4);
      char c;
      if(randChar == 0) c='A';
      else if(randChar == 1) c='C';
      else if(randChar == 2) c= 'G';
      else c='T';
      read[randInd] = c;
    }
    file << read << endl;
  }
  file.close();

  ofstream file2;
  file2.open("rands_debug.txt");
  for(auto read : rands){
    file2 << read << endl;
  }
  file2.close();

  ofstream file3;
  file2.open("myseq_debug.txt");
  file2 << string(myseq) << endl;
  file2.close();
  return rtn;
}

void BruteForceRecon(string ref, vector<string> &shortReads){
  // reconstruct DNA
  // init
  clock_t start = clock();
  int k = shortReads[0].length();
  int n = shortReads.size();
  char recon[ref.length()];
  cout << "init done!" << endl;
  for(int i=0; i<ref.length(); i++) recon[i]='A';
  // matching
  for(int readInd =0; readInd<shortReads.size(); readInd++){
    if(readInd % 1000 == 0) cout << readInd << " progress!" << endl;
    string read = shortReads[readInd];
    for(int j=0; j<ref.length()-k-1; j++){
      //sliding window
      int diff = 0;
      // diff calculation
      for(int ind=0; ind<k; ind++){
        if(diff > 2) break;
        if(read[ind] != ref[j+ind]){
          diff += 1;
        }
      }
      // match
      if(diff <= 2){
        for(int i=0; i<k; i++) recon[j+i] = read[i];
      }
    }
  }
  // matching accuracy
  int matched = 0;
  for(int i=0; i<ref.length(); i++){
    if(ref[i] == recon[i]) matched += 1;
  }
  clock_t end = clock();
  double duration = (double)(end - start) / CLOCKS_PER_SEC;
  cout << " Time spend        :: " << duration << endl;
  cout << " matching accuracy :: " << matched << "/" << ref.length() << endl;
  return ;
}

int overlap(string s1, string s2,int min_length){
  int start = 0;
  //cout << s1 << " " << s2 << endl;
  start = s1.find(s2.substr(0,min_length));
  //cout << start << endl;
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
    }
    return start;
  }
  return start;
}

vector <int> pick_maximal_overlap(vector <string> reads, int k){
  vector <int> out;
  int best_olen = 0;
  int reada, readb;

  //cout << " pick maximal overlap :: reads " << endl;
  //for(auto s : reads)cout << " " << s << endl;
  //cout << "=====" << endl;

  map <string, vector<string> > index;
  for(auto read : reads){
    for(int i=0; i<read.length()-k+1; i++){
      index[read.substr(i,k)].push_back(read);
      //cout << read.substr(i,k) << " ";
    }//cout << endl;
  }

  int i,j;i=0;j=0;
  for(auto r : reads){
    //cout << r << " " << r.substr(r.length()-k,k) << endl;
    for(auto o : index[r.substr(r.length()-k,k)] ){
      if(r != o){
        //cout << " overlap find " << r << " " << o << endl;
        int olen = overlap(r,o,o.length());
        if(olen > best_olen){
          vector <string>::iterator itera = find(reads.begin(), reads.end(), r);
          vector <string>::iterator iterb = find(reads.begin(), reads.end(), o);
          reada = distance(reads.begin(), itera);
          readb = distance(reads.begin(), iterb);
          //cout << "best olen " << olen<< " "<< r << " "<<o<< " " << reada << " " << readb << endl;
          best_olen = olen;
        }
      }
    }
  }
  out.push_back(reada);
  out.push_back(readb);
  out.push_back(best_olen);
  return out;
}

string greedySCS(vector <string> shortReads, int k){
  clock_t start = clock();
  vector<int> out = pick_maximal_overlap(shortReads, k);
  int reada = out[0]; int readb = out[1]; int olen = out[2];
  int count = 0;
  while(olen>0){
    cout << count << " "<< shortReads.size() << endl;
    //cout << reada << " "<< readb << " "<< olen << endl;
    string read_a = shortReads[reada];
    string read_b = shortReads[readb];
    //cout << read_a << " " << read_b << endl;
    if(reada > readb){
      shortReads.erase(shortReads.begin()+reada);
      shortReads.erase(shortReads.begin()+readb);
    }else{
      shortReads.erase(shortReads.begin()+readb);
      shortReads.erase(shortReads.begin()+reada);
    }
    shortReads.push_back(read_a.substr(0, olen) + read_b);
    vector<int> tmp = pick_maximal_overlap(shortReads, k);
    reada = tmp[0]; readb = tmp[1]; olen = tmp[2];
    count += 1;
  }
  string ret = "";
  for(auto s : shortReads){
    ret += s;
  }
  clock_t end = clock();
  cout << "Time spend :: "<< end-start<<endl;
  return ret;
}

string ShortestCommonSuperstring(vector<string> &shortReads){
  bool print = true;
  string shortest_sup = "x";
  int n = shortReads.size();
  vector<int> perm;
  for(int i=0; i<n; i++) perm.push_back(i);
  sort(perm.begin(), perm.end());
  clock_t start = clock();
  int count =0;
  do {
      cout << count << endl;
      count++;
      string sup = shortReads[perm[0]];
      if(print) cout << "sup:: "<< sup << " "<< perm[0]<<perm[1] << perm[2]<<endl;
      if(print) cout << sup << endl;
      int indent =0;
      for(int i=0; i<n-1; i++){
        string s1 = shortReads[perm[i]];
        string s2 = shortReads[perm[i+1]];
        int olen = overlap(s1, s2, s2.length());
        if(olen==-1){
          indent += s1.length();
          sup += s2;
        }
        else{
          indent += olen;
          sup += s2.substr(s1.length()-olen,s2.length()-s1.length()+olen);
        }
        if(print) for(int i=0; i<indent; i++) cout << " ";
        if(print) cout << s2 << endl;
      }
      if(print) cout << sup << " "<< sup.length() << endl;
      if(shortest_sup =="x" || sup.length() < shortest_sup.length()){
        shortest_sup = sup;
         if(print) cout << "shortest :: "<< shortest_sup << " "<<sup.length()<< endl;
      }
      if(print) cout << endl;
    }while (next_permutation(perm.begin(), perm.end()));
  clock_t end = clock();
  cout << "Time spend :: "<< end-start<<endl;
  return shortest_sup;
}

void denovo(vector<string> &shortReads){
  cout << "denovo" << endl;
  map <string, vector<string> > G;
  int k = shortReads[0].length();
  for(auto read : shortReads){
    // read :: "ATC...TGC"
    string left = read.substr(0,k-1);
    string right = read.substr(1,k-1);

    G[left].push_back(right);

  }

  // Graph Print
  for(auto read : shortReads){
    string left = read.substr(0,k-1);
    string right = read.substr(1,k-1);

    cout << left << " -> ";
    for(int i=0; i<G[left].size(); i++){
      cout << G[left][i] << " ";
    }
    cout << endl;
  }


  return;
}

int main(){
  init(10000); // random Sequence generate
  reviseSeq(); // duplicate removal

  string ref = import();
  vector<string> shortReads = Sequencer(400, 30, ref);
  //string arr[] = {"ACGGATGAGC", "GAGCGGA", "GAGCGAG"};
  // Initialize vector with a string array
  //vector<string> shortReads(arr, arr + sizeof(arr)/sizeof(string));
  //BruteForceRecon(ref, shortReads);
  string recon = greedySCS(shortReads, 2);
  cout << recon.length() << endl;
  int matched = 0;
  for(int i=0; i<ref.length(); i++){
    if(ref[i] == recon[i]) matched += 1;
  }
  cout << matched << "/" << ref.length() << endl;
  return 0;
}
