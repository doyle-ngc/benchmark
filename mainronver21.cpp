#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <queue>
#include <stack>
#include <unistd.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <tuple>
#include <cmath>
#include <time.h> 
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/foreach.hpp>
#include <cstdio>

using namespace boost;
using namespace std;


struct Cnot{
  int cbit_no,tbit_no;
};

struct Node{
  int id,x,y,k;
};

class state{
  public:
  int w;
  int h;
  vector<vector <int> > grid;
  int depth=0;
  float evaluation;
  bool operator<(const state &s) const {
    return evaluation > s.evaluation;
  }
  state(int w0, int h0){
    w = w0;
    h = h0;
    grid.resize(w*h);
    for(int i=0;i<w;i++){
      grid[i].resize(h);
    }
  }
};

typedef adjacency_list<listS, vecS, directedS,Cnot> Graph;
typedef pair<int, int> Edge;
typedef graph_traits<Graph>::vertex_descriptor   Vertex;

int w,h,qubit_no;
vector<int> evaluate_progress;


void inputFileGenerate(int gate_number, int grid_size){
  fstream outputfile("test.txt",ios::out);
   if(outputfile.fail()) cout<<"file open faild"<<endl;
   int check_gate_number = 0;
    srand((unsigned) time(NULL));
    while(1){
      if(check_gate_number == gate_number) break;
      else{
        int x = rand()%grid_size;
        int y = rand()%grid_size;
        if(x != y){
          outputfile<<y<<","<<x<<std::endl;
          cout<<y<<","<<x<<endl;
          check_gate_number++;
        }
      }
    }
    outputfile.close();
}

void benchmarkFileGenerate(vector<state> states,int file_number){

char fname[30];
sprintf(fname, "100qubits_layouts_%d.txt",file_number);
ofstream file(fname);

  fstream outputfile(fname,ios::out);
  if(outputfile.fail()) cout<<"outputfile open faild"<<endl;

  outputfile<<w<<" "<<h<<" "<<states.size()<<endl;;
  for(auto e:states){
    for(int j=0;j<h;j++){
      for(int i=0;i<w;i++){
        outputfile<<j<<" "<<i<<" "<<e.grid[i][j]<<endl;;
      }
    }
  }

  outputfile.close();

}

int circuitGenerate(vector< pair<Cnot,int> > &circuit,string input){
  string str;
  int qubit_no = 0;
  int gate_no = 1;
  fstream ifs(input,std::ios::in);
    if(ifs.fail()) cout<<"file open failed txt"<<endl;
    while(getline(ifs, str)){
      int cbit = 0;
      int tbit = 0;
      sscanf(str.data(), "%d,%d", &cbit, &tbit);
      pair<Cnot,int> gate;
      Cnot cnot;
      cnot.cbit_no = cbit;
      cnot.tbit_no = tbit;
      if(cbit>qubit_no){
        qubit_no = cbit;
      }
      if(tbit>qubit_no){
        qubit_no = tbit;
      }
      gate.first = cnot;
      gate.second = gate_no;
      circuit.push_back(gate);
      gate_no++;
    }

    ifs.close();

    return (qubit_no+1);
}

void nodeGenerate(int x,int y,int k,vector<Node> &node_list){

  //node_list.clear();

  int node_no = 1;
  for(int i=0;i<k;i++){
    for(int j=0;j<y;j++){
      for(int l=0;l<x;l++){
        Node node; 
        node.id = node_no;
        node.x = l;
        node.y = j;
        node.k = i;
        node_list.push_back(node);
        node_no++;
      }
    }
  }

}

int getNodeNo(int x,int y,int k,vector<Node> node_list){
  for(auto e:node_list){
    if(e.x == x && e.y == y && e.k == k){
      return e.id;
    }
  }
  return -1;
}

int getXNo(int node_no,vector<Node> node_list){
  for(auto e:node_list){
    if(e.id == node_no){
      return e.x;
    }
  }
  return -1;
}

int getYNo(int node_no,vector<Node> node_list){
  for(auto e:node_list){
    if(e.id == node_no){
      return e.y;
    }
  }
  return -1;
}

int getKNo(int node_no,vector<Node> node_list){
  for(auto e:node_list){
    if(e.id == node_no){
      return e.k;
    }
  }
  return -1;
}

void addCnfFixedPair(int cons1,int cons2, vector< string > &cnf_fixed){
  string str;
  str += ("-");
  str += std::to_string(cons1);
  str += (" ");
  str += ("-");
  str += std::to_string(cons2);
  str += (" ");
  str += ("0");
  cnf_fixed.push_back(str);
}

void addCnfOnlySum(vector< string > &cnf_fixed,int &cnf_fixed_no){

  for(int i=0;i<qubit_no;i++){
    string str;
    for(int j=1;j<=w*h;j++){   
      str += to_string((w*h*i)+j);
      str += (" ");
    }
    str += ("0");
    cnf_fixed_no++;
    cnf_fixed.push_back(str);
  }

}

void addCnfOnlyPair(vector< Node > node_list,vector< string > &cnf_fixed,int &cnf_fixed_no){
  
  for(int k=0;k<qubit_no;k++){
    for(int j=0;j<h;j++){
      for(int i=0;i<w;i++){
        for(int m=0;m<h;m++){
          for(int l=0;l<w;l++){
            if(getNodeNo(i,j,k,node_list)<getNodeNo(l,m,k,node_list)){
              addCnfFixedPair(getNodeNo(i,j,k,node_list),getNodeNo(l,m,k,node_list),cnf_fixed);
              cnf_fixed_no++;
            }
          }
        }
      }
    }
  }

}

void addCnfAtmostPair(vector< string > &cnf_fixed,vector< Node > node_list,int &cnf_fixed_no){

  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      for(int k=0;k<qubit_no;k++){
        for(int l=k+1;l<qubit_no;l++){
          addCnfFixedPair(getNodeNo(i,j,k,node_list),getNodeNo(i,j,l,node_list),cnf_fixed);
          cnf_fixed_no++;
        }
      }
    }
  }

}

void addCnfAdjacent(int c,int t,vector< Node > node_list,vector< string > &cnf_adjacent,int &cnf_adjacent_no){

  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      string str;
      str += ("-");
      str += to_string(getNodeNo(i,j,c,node_list));
      str += (" ");
      if(i-1!=-1) {
        str += to_string(getNodeNo(i-1,j,t,node_list));
        str += (" ");
      }
      if(i+1!=w) {
        str += to_string(getNodeNo(i+1,j,t,node_list));
        str += (" ");
      }
      if(j-1!=-1) {
        str += to_string(getNodeNo(i,j-1,t,node_list));
        str += (" ");
      }
      if(j+1!=h) {
        str += to_string(getNodeNo(i,j+1,t,node_list));
        str += (" ");
      }
      str += ("0");
      cnf_adjacent_no++;
      cnf_adjacent.push_back(str);
    }
  }

}

void addCnfSearchedAnswer(vector<int> result_int,vector< string > &cnf_searched_answer,int &cnf_searched_answer_no){
  string str;
  for(auto e:result_int){
    str += ("-");
    str += to_string(e);
    str += (" ");
  }
  str += ("0");
  cnf_searched_answer.push_back(str);
  cnf_searched_answer_no++;
}

void addCnfVacancySum(vector< string > &cnf_fixed,int &cnf_fixed_no,vector<int> vacancy_grid){

  for(auto e:vacancy_grid){   
    for(int i=0;i<qubit_no;i++){
      string str;
      str += ("-");
      str += to_string(e+(i*w*h));
      str += (" ");
      str += ("0");
      cnf_fixed_no++;
      cnf_fixed.push_back(str);
    } 
  }

}

void deleteCnfAdjacent(vector< string > &cnf_adjacent,int &cnf_adjacent_no){

  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      cnf_adjacent.pop_back();
      cnf_adjacent_no--;
    }
  }

}

void fixedCnfEncoder(vector<Node> node_list,vector<string> &cnf_fixed,int &cnf_fixed_no,vector<int> vacancy_grid){

  cnf_fixed_no = 0;
  cnf_fixed.clear();

  addCnfOnlySum(cnf_fixed,cnf_fixed_no);
  addCnfOnlyPair(node_list,cnf_fixed,cnf_fixed_no);
  addCnfAtmostPair(cnf_fixed,node_list,cnf_fixed_no);
  addCnfVacancySum(cnf_fixed,cnf_fixed_no,vacancy_grid);

}

void CnfstartEncoder(vector<string> &cnf_start,int total_cnf_no){
  cnf_start.clear();
  //total_cnf_no = cnf_fixed_no + cnf_adjacent_no;
  string str;
  str += ("p");
  str += (" ");
  str += ("cnf");
  str += (" ");
  str += to_string(qubit_no*w*h);
  str += (" ");
  str += to_string(total_cnf_no);
  cnf_start.push_back(str);
}

void CnfInputGenerate(vector<string> &cnf_start,vector<string> cnf_fixed,vector<string> cnf_adjacent,vector<string> cnf_searched_answer,int total_cnf_no){

  CnfstartEncoder(cnf_start,total_cnf_no);
  fstream outputfile("input.dimacs",ios::out);
  if(outputfile.fail()) cout<<"file open faild"<<endl;

  for(auto e:cnf_start){
    outputfile<<e<<endl;
  }

  for(auto e:cnf_fixed){
    outputfile<<e<<endl;
  }

  for(auto e:cnf_adjacent){
    outputfile<<e<<endl;
  }

  for(auto e:cnf_searched_answer){
    outputfile<<e<<endl;
  }

  outputfile.close();
}

void SAT(vector<string> &cnf_start,vector<string> &cnf_fixed,vector<string> &cnf_adjacent,vector<string> &cnf_searched_answer,int &total_cnf_no,int cnf_fixed_no,int cnf_adjacent_no,int cnf_searched_answer_no){

  int pid,code,status;
  pid_t result;

  total_cnf_no = cnf_fixed_no+cnf_adjacent_no+cnf_searched_answer_no;

  CnfInputGenerate(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no);
  
  pid = fork();
      //std::cout<<"pid="<<pid<<std::endl;
      if(pid == -1){
        cout<<"pid error"<<endl;
      }
      if(pid == 0){
        execl("./glueminisat-simp","./glueminisat-simp","input.dimacs","output",NULL);
      }
      else{
        result = wait(&status);
      }

      if(result < 0){
        fprintf(stderr, "Wait Error.\n\n");
        exit(-1);
      }
      
      // 終了ステータスのチェック
      if(WIFEXITED(status)){
      //printf("子プロセス終了");
      code = WEXITSTATUS(status);
      //printf("コード : %d\n", code);
      }
      else{
        cout<<"wait失敗"<<endl;
        printf("終了コード : %d\n", status);
      }

}

void getCnfAnswer(vector<int> &result_int){

  result_int.clear();

  fstream ifs("output",ios::in);
  if(ifs.fail()) cout<<"file open failed"<<endl;
  string str;
  int temp;
  
  while(!ifs.eof()){
    ifs>>str;
    if(str!="SAT"){
      temp = stoi(str);
      if(temp>0){
        result_int.push_back(temp);
      }
    }
  }

  ifs.close();

}

int checkSAT(){
  fstream ifs("output",ios::in);
  if(ifs.fail()) cout<<"file open failed"<<endl;
  string str;

  while(!ifs.eof()){
    ifs>>str;
    if(str=="SAT"){
      return 1;
    }
    else if(str=="UNSAT"){
      return -1;
    }  
  }

  return 0;
}

bool isSameState(state s,state goal){
  for(int j=0;j<h;j++){
    for(int i = 0;i<w;i++){
      if(s.grid[i][j] != goal.grid[i][j]){
        return false;
      }
    }
  }
  return true;
}

bool isSameCircuit(vector<Cnot> circuit1, vector<Cnot> circuit2){

    if(circuit1.size()!=circuit2.size()) return false;
    else{
        auto itr2 = circuit2.begin();
        for(auto itr1 = circuit1.begin(); itr1!=circuit1.end(); itr1++){
            auto gate1 = *itr1;
            auto gate2 = *itr2;
            if(gate1.cbit_no!=gate2.cbit_no||gate1.tbit_no!=gate2.tbit_no) return false;
            itr2++;
        }
    }

    return true;
}

bool checkAdjacentCondition(state s, vector<Cnot> subcircuit){

  vector<int> check_adjacents;
  int result = 1;

  for(auto e:subcircuit){
    int check_adjacent = 0;
    for(int j=0;j<h;j++){
      for(int i=0;i<w-1;i++){
        if((s.grid[i][j]==e.cbit_no&&s.grid[i+1][j]==e.tbit_no)||(s.grid[i][j]==e.tbit_no&&s.grid[i+1][j]==e.cbit_no)){
          check_adjacent = 1;
        }
      }
    }
    for(int j=0;j<h-1;j++){
      for(int i=0;i<w;i++){
        if((s.grid[i][j]==e.cbit_no&&s.grid[i][j+1]==e.tbit_no)||(s.grid[i][j]==e.tbit_no&&s.grid[i][j+1]==e.cbit_no)){
          check_adjacent = 1;
        }
      }
    }
    check_adjacents.push_back(check_adjacent);
  }
  for(auto e:check_adjacents){
    if(e==0){
      result = 0;
      break;
    }
  }

  return result;
}

int manhattanDist(state s,state goal){

  //int temp_max_dist = 0;
  int eval = 0;
  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      for(int l=0;l<h;l++){
        for(int k=0;k<w;k++){
          if(s.grid[i][j] == goal.grid[k][l] && s.grid[i][j] < qubit_no){
            eval += (abs(i-k)+abs(j-l));
          }
        }
      }
    }
  }
  return eval;
}

int differentNodeNo(state s,state goal){
  int different_node_no = 0;
  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      if(s.grid[i][j] != goal.grid[i][j] && s.grid[i][j] < qubit_no){
        different_node_no++;
      }
    }
  }
  return different_node_no;
}

int manhattanDist2(state s,vector<Cnot> subcircuit){

  int eval = 0;
  for(auto e:subcircuit){
    for(int i=0;i<w;i++){
      for(int j=0;j<h;j++){
        for(int l=0;l<h;l++){
          for(int k=0;k<w;k++){
            if((e.cbit_no==s.grid[i][j]&&e.tbit_no==s.grid[k][l])||(e.tbit_no==s.grid[i][j]&&e.cbit_no==s.grid[k][l])){
              eval += (abs(i-k)+abs(j-l));
            }
          }
        }
      }
    }
  }

  return eval;

}

void printState(state s){
  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      cout<<s.grid[i][j]<<" ";
    }
    cout<<endl;
  }
  cout<<endl;
}

void deleteSubcircuitgates(vector<pair<Cnot,int>> &circuit,vector<Cnot> &subcircuit){

  for(auto e:subcircuit){
    for(auto gate_itr = circuit.begin();gate_itr!=circuit.end();gate_itr++){
       auto gate = *gate_itr;
       int check = 0;
       if(e.cbit_no==gate.first.cbit_no&&e.tbit_no==gate.first.tbit_no&&check==0){
         //cout<<"erase"<<gate.first.cbit_no<<","<<gate.first.tbit_no<<","<<gate.second<<endl;
         circuit.erase(gate_itr);
         gate_itr--;
         check=1;
    }
    }
  }
}


int AStar_re(state goal,priority_queue<state> &open,vector<state> &closed){
  int check_searched = 0;
  int swap_no;
  if(open.empty()){
    cout<<"open list error"<<endl;
  }
  state now = open.top();
  open.pop();
  /*std::cout<<"now state"<<std::endl;
  printState(now);
  std::cout<<"open list size:"<<open.size()<<std::endl;
  std::cout<<"closed list size:"<<closed.size()<<std::endl;
  std::cout<<"searched list size:"<<searched_state.size()<<std::endl;
  std::cout<<std::endl;*/
  if(isSameState(now,goal)){
    cout<<"finished"<<endl;
    printState(now);
    cout<<"swap number:"<<now.depth<<endl;
    swap_no = now.depth;
  }  
  else{
    state s = now;
    for(auto e:closed){
      if(isSameState(e,s) == true){
        check_searched = 1;
        //<<"check_searched"<<check_searched<<std::endl;
        }
    }
    if(check_searched == 0){
      for(int j=0;j<h;j++){
        for(int i=0;i<w-1;i++){
          if(s.grid[i][j]<qubit_no&&s.grid[i+1][j]<qubit_no){
            swap(s.grid[i][j],s.grid[i+1][j]);
            printState(s);
            s.depth++;
            s.evaluation =  s.depth + manhattanDist(s,goal);// + differentNodeNo(s,goal)  ;
            /*std::cout<<"dist:"<<manhattanDist(s,goal)<<std::endl;
            <<"depth:"<<s.depth<<std::endl;
            std::cout<<"evaluation:"<<s.evaluation<<std::endl;
            std::cout<<std::endl;*/
            //if(open.size()<4000000){
              open.push(s);
            //}
            s = now;
          }
        }
      }
      for(int i=0;i<w;i++){
        for(int j=0;j<h-1;j++){
          if(s.grid[i][j]<qubit_no&&s.grid[i][j+1]<qubit_no){
          swap(s.grid[i][j],s.grid[i][j+1]);
          //printState(s);
          s.depth++;
          s.evaluation = s.depth + manhattanDist(s,goal);//+ differentNodeNo(s,goal);
          /*std::cout<<"dist:"<<manhattanDist(s,goal)<<std::endl;
          std::cout<<"depth:"<<s.depth<<std::endl;
          std::cout<<"evaluation:"<<s.evaluation<<std::endl;
          std::cout<<std::endl;*/
          //if(open.size()<4000000){
            open.push(s);
          //}
          s = now;
          }
        }
      }
      closed.push_back(s);
    }
    swap_no=AStar_re(goal,open,closed);
  }

  return(swap_no);
}

int Astar2_re(state &start, std::vector<Cnot> subcircuit,priority_queue<state> &open,vector<state> &closed){

  int check_searched = 0;
  int swap_no;

  if(open.empty()){
    cout<<"open list error"<<endl;
  }

  state now = open.top();
  open.pop();
  if(checkAdjacentCondition(now,subcircuit)){
    cout<<"finished"<<endl;
    cout<<"swap number:"<<now.depth<<endl;
    swap_no = now.depth;
    start = now;
  }
  else{
    state s = now;
    for(auto e:closed){
      if(isSameState(s,e)== true){
        check_searched = 1;
      }
    }
    if(check_searched == 0){
      for(int j=0;j<h;j++){
        for(int i=0;i<w-1;i++){
          if(s.grid[i][j]<qubit_no&&s.grid[i+1][j]<qubit_no){
            swap(s.grid[i][j],s.grid[i+1][j]);
            s.depth++;
            s.evaluation = manhattanDist2(s,subcircuit)  + s.depth;
            cout<<"dist"<<manhattanDist2(s,subcircuit)<<endl;
            cout<<"depth:"<<s.depth<<endl;
            cout<<"evaluation:"<<s.evaluation<<endl;
            cout<<endl;
            open.push(s);
            s = now;
          }
        }
      }
      for(int i=0;i<w;i++){
        for(int j=0;j<h-1;j++){
          if(s.grid[i][j]<qubit_no&&s.grid[i][j+1]<qubit_no){
            std::swap(s.grid[i][j],s.grid[i][j+1]);
            s.depth++;
            s.evaluation = manhattanDist2(s,subcircuit) + s.depth;
            open.push(s);
            s = now;
          }
        }
      }
      closed.push_back(s);
    }
    swap_no = Astar2_re(start,subcircuit,open,closed);
  }

  return swap_no;

}

int Astar(state start, state goal){

  int swap_no;
  priority_queue<state> open;
  vector<state> closed;
  //std::vector<state> searched_state;

  start.evaluation = manhattanDist(start,goal);//+differentNodeNo(start,goal);
  open.push(start);
  swap_no = AStar_re(goal,open,closed);

  return(swap_no);
  
}

int Astar2(state &start, vector<Cnot> subcircuit){

  int swap_no;
  priority_queue<state> open;
  vector<state> closed;

  start.evaluation = manhattanDist2(start, subcircuit);
  open.push(start);
  swap_no = Astar2_re(start,subcircuit,open,closed);

  return(swap_no);

}

void generateState(state &s,vector<Node> node_list,vector<int> result_int,vector<int> vacancy_grid){
  

  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      s.grid[i][j] = qubit_no;
    }
  }

  for(auto e:result_int){
    s.grid[getXNo(e,node_list)][getYNo(e,node_list)] = getKNo(e,node_list);
  }
}

void adjacentPairs(state s,vector< pair<int,int> > &adjacent_pairs,vector<Node> node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,vector<Cnot> subcircuit){

  adjacent_pairs.clear();
  vector<int> moves;

  for(int j=0;j<h;j++){
    for(int i=0;i<w-1;i++){
      if(s.grid[i][j]<qubit_no && s.grid[i+1][j]<qubit_no){
        pair<int,int>  adjacent_pair;
        adjacent_pair.first = s.grid[i][j];
        adjacent_pair.second = s.grid[i+1][j];
        adjacent_pairs.push_back(adjacent_pair);
      }
    }
  }

  for(int j=0;j<h-1;j++){
    for(int i=0;i<w;i++){
      if(s.grid[i][j]<qubit_no && s.grid[i][j+1]<qubit_no){
        pair<int,int>  adjacent_pair;
        adjacent_pair.first = s.grid[i][j];
        adjacent_pair.second = s.grid[i][j+1];
        adjacent_pairs.push_back(adjacent_pair);
      }
    }
  }

  for(auto e:subcircuit){
    int check_cbit=0;
    int check_tbit=0;
    for(auto f:moves){
      if(f==e.cbit_no) check_cbit=1;
      if(f==e.tbit_no) check_tbit=1;
    }
    if(check_cbit==0) moves.push_back(e.cbit_no);
    if(check_tbit==0) moves.push_back(e.tbit_no);
  }

  for(int j=0;j<h;j++){
    for(int i=0;i<w;i++){
      for(auto e:subcircuit){
        if(s.grid[i][j]==e.cbit_no||s.grid[i][j]==e.tbit_no){
          if(i-1!=-1) moves.push_back(s.grid[i-1][j]);
          if(i+1!=w) moves.push_back(s.grid[i+1][i]);
          if(j-1!=-1) moves.push_back(s.grid[i][j-1]);
          if(j+1!=h)  moves.push_back(s.grid[i][j+1]);
        }
      }
    }
  }

  for(auto e:adjacent_pairs){
        int check_adjacent = 0;
        for(auto j:moves){
          if(j==e.first || j==e.second){
            check_adjacent = 1;
          }
        }
        
        if(check_adjacent==0){
          addCnfAdjacent(e.first,e.second,node_list,cnf_adjacent,cnf_adjacent_no);
        }
      }
}

void searchBestPlacement(vector<int> &result_int,vector<string> &cnf_start,vector<string> &cnf_fixed, vector<string> cnf_adjacent, int cnf_fixed_no,int cnf_adjacent_no, vector<int> vacancy_grid, state start, state &goal,vector<Node> node_list){

  vector <string> cnf_searched_answer;
  vector<int> temp_result_int;
  int cnf_searched_answer_no = 0;
  int SAT_count = 0;
  int state_eval,temp_state_eval,total_cnf_no;
  state temp_goal(w,h);

  temp_state_eval = 1000000;


  while(checkSAT()==1&&SAT_count<500){

        getCnfAnswer(result_int);
        addCnfSearchedAnswer(result_int,cnf_searched_answer,cnf_searched_answer_no); 
        generateState(temp_goal,node_list,result_int,vacancy_grid);
        state_eval = manhattanDist(start,temp_goal) ;//+ differentNodeNo(start,temp_goal);
        //manhattan_dists.push_back(state_eval);
        if(state_eval<temp_state_eval){
          temp_state_eval = state_eval;
          goal = temp_goal;
          temp_result_int = result_int;
          evaluate_progress.push_back(temp_state_eval);
        }
        SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
        SAT_count++;
        cout<<SAT_count<<endl;
    }


}

void generateDependencyGraph(Graph &graph,std::vector< pair<Cnot,int> > circuit,vector< vector <int> > connect){

   graph.clear();
  int gate_count = 0;
  map<int, Graph::vertex_descriptor> desc;
  desc[0] = add_vertex(graph);
  for(auto e:circuit){
       desc[e.second] = add_vertex(graph);
       graph[desc[e.second]].cbit_no = e.first.cbit_no;
       graph[desc[e.second]].tbit_no = e.first.tbit_no;
  }
  
  for(auto outer_itr = circuit.begin();outer_itr!=circuit.end();outer_itr++){  //generate a dependency graph
    auto outer_gate = *outer_itr;
    //std::cout<<outer_gate.first.cbit_no<<" "<<outer_gate.first.tbit_no<<" "<<outer_gate.second<<std::endl;
    for(auto inner_itr = outer_itr+1;inner_itr!=circuit.end();inner_itr++){
      auto inner_gate = *inner_itr;
      if((outer_gate.first.cbit_no==inner_gate.first.tbit_no)||(outer_gate.first.tbit_no==inner_gate.first.cbit_no)){
        //std::cout<<outer_gate.first.cbit_no<<" "<<outer_gate.first.tbit_no<<" "<<outer_gate.second<<std::endl;
        //std::cout<<inner_gate.first.cbit_no<<" "<<inner_gate.first.tbit_no<<" "<<inner_gate.second<<std::endl;
        //std::cout<<std::endl;
          add_edge(desc[outer_gate.second],desc[inner_gate.second],graph);
          connect[outer_gate.second][inner_gate.second] = 1; 
      }
    }
    gate_count++;
  }

  for(int j=1;j<=gate_count;j++){  //add edges from node0
    int check_edge = 0;
    for(int i=1;i<=gate_count;i++){
      if(connect[i][j] == 1) check_edge = 1;
    } 
    if(check_edge == 0){
      add_edge(desc[0],desc[j],graph);
      connect[0][j] = 1;
    }
  }
  
  for(int i=1;i<=gate_count;i++){  //remove reduntant edges
    for(int j=1;j<=gate_count;j++){
      if(connect[i][j]==1) {
        std::vector<boost::default_color_type> color(num_vertices(graph), boost::white_color);
        if(is_reachable(desc[i], desc[j], graph, color.data())==0) {
          add_edge(desc[i],desc[j],graph);
        }
      }
    }
  }

  auto vertices = adjacent_vertices(desc[0],graph);
  cout<<"foreach"<<endl;
  BOOST_FOREACH(auto e,vertices){
    cout<<e<<endl;
  }

}

void generateConstructedDependencyGraph(Graph &graph,std::vector<std::pair<Cnot,int>> &circuit,vector<std::vector <int>> &connect,int gate_no,std::vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,std::vector<Cnot> &subcircuit){

  graph.clear();
  subcircuit.clear();
  cnf_adjacent.clear();
  map<int, Graph::vertex_descriptor> desc;
  desc[0] = add_vertex(graph);
  vector<int> removed_gates_id;
  for(auto e:circuit){
       desc[e.second] = add_vertex(graph);
       graph[desc[e.second]].cbit_no = e.first.cbit_no;
       graph[desc[e.second]].tbit_no = e.first.tbit_no;
  }
  int added_gate_count = 0;
  vector<int> added_vertices;
  vector<Graph> graphs;
  
  for(auto outer_itr = circuit.begin();outer_itr!=circuit.end();++outer_itr){  //generate a dependency graph
    auto outer_gate = *outer_itr;
    //std::cout<<outer_gate.first.cbit_no<<" "<<outer_gate.first.tbit_no<<" "<<outer_gate.second<<std::endl;
    for(auto inner_itr = outer_itr+1;inner_itr!=circuit.end();++inner_itr){
      auto inner_gate = *inner_itr;
      if((outer_gate.first.cbit_no==inner_gate.first.tbit_no)||(outer_gate.first.tbit_no==inner_gate.first.cbit_no)){
          int check_added_outer,check_added_inner;
          check_added_outer=check_added_inner=0;
          auto vertices_graph = vertices(graph);
          add_edge(desc[outer_gate.second],desc[inner_gate.second],graph);
          connect[outer_gate.second][inner_gate.second] = 1; 
          removed_gates_id.push_back(desc[outer_gate.second]);
          removed_gates_id.push_back(desc[inner_gate.second]);
          for(auto e:added_vertices){
            if(e==outer_gate.second) check_added_outer = 1;
            if(e==inner_gate.second) check_added_inner = 1;
          }
          if(check_added_outer==0) added_vertices.push_back(outer_gate.second);
          if(check_added_inner==0) added_vertices.push_back(inner_gate.second);
          added_gate_count++;
          //cout<<outer_gate.second<<inner_gate.second<<endl;
      }
      if(added_vertices.size()>=gate_no) break;
    }
    if(added_vertices.size()>=gate_no) break;
  }

  /*for(int i=1;i<=num_vertices(graph);i++){
    int check_vertices = 0;
    for(int j=1;j<=num_vertices(graph);j++){
      if(connect[i][j]==1||connect[j][i]==1) check_vertices = 1;
    }
    if(check_vertices==0) remove_vertex(desc[i],graph);
  }*/

  graphs.push_back(graph);

  /*for(int j=1;j<num_vertices(graph);j++){  //add edges from node0
    int check_edge = 0;
    for(int i=1;i<num_vertices(graph);i++){
      if(connect[i][j]) check_edge = 1;
    } 
    if(check_edge == 0){
      add_edge(desc[0],desc[j],graph);
      connect[0][j] = 1;
      added_gate_count++;
    }
  }
  
  for(int i=1;i<added_gate_count;i++){  //remove reduntant edges
    for(int j=1;j<added_gate_count;j++){
      if(connect[i][j]==1){
        remove_edge(desc[i],desc[j],graph);
        std::vector<boost::default_color_type> color(num_vertices(graph), boost::white_color);
        if(is_reachable(desc[i], desc[j], graph, color.data())==0){
          add_edge(desc[i],desc[j],graph);
        }
      }
    }
  }*/

  for(auto e:added_vertices){
    Cnot cnot;
    cnot.cbit_no = graph[desc[e]].cbit_no;
    cnot.tbit_no = graph[desc[e]].tbit_no;
    subcircuit.push_back(cnot);
    addCnfAdjacent(graph[desc[e]].cbit_no,graph[desc[e]].tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }

  cout<<"added_vertices.size"<<added_vertices.size()<<endl;

  if(added_vertices.size()==0){
    for(auto e:circuit){
    Cnot cnot;
    cnot.cbit_no = e.first.cbit_no;
    cnot.tbit_no = e.first.tbit_no;
    subcircuit.push_back(cnot);
    addCnfAdjacent(e.first.cbit_no,e.first.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }
  }

  
  for(auto gate_itr = circuit.begin();gate_itr!=circuit.end();gate_itr++){
    auto gate = *gate_itr;
    int check = 0;
    for(auto e:removed_gates_id){
      if(e==gate.second&&check==0) {
        circuit.erase(gate_itr);
        gate_itr--;
        check=1;
      }
    }
  }

}

void generateConstructedDependencyGraph2(Graph &graph,std::vector<std::pair<Cnot,int>> &circuit,vector<std::vector <int>> &connect,int gate_no,std::vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,std::vector<Cnot> &subcircuit){

  graph.clear();
  subcircuit.clear();
  cnf_adjacent.clear();
  map<int, Graph::vertex_descriptor> desc;
  desc[0] = add_vertex(graph);
  vector<int> removed_gates_id;
  for(auto e:circuit){
       desc[e.second] = add_vertex(graph);
       graph[desc[e.second]].cbit_no = e.first.cbit_no;
       graph[desc[e.second]].tbit_no = e.first.tbit_no;
  }
  int added_gate_count = 0;
  vector<int> added_vertices;
  vector<Graph> graphs;
  vector< pair<Cnot,int> > not_added_gates,constructed_subcircuit;

  for(auto gate_itr = circuit.begin();gate_itr!=circuit.end();++gate_itr){

    auto gate = *gate_itr;
    int check_add = 0;
    for(auto e:not_added_gates){
      if(gate.first.cbit_no==e.first.tbit_no||gate.first.tbit_no==e.first.cbit_no) check_add = 1;
    }
    if(check_add==0) constructed_subcircuit.push_back(gate);
    else not_added_gates.push_back(gate);
    if(constructed_subcircuit.size()>=gate_no) break;
  }

  for(auto e:constructed_subcircuit){
    Cnot cnot;
    cnot.cbit_no = e.first.cbit_no;
    cnot.tbit_no = e.first.tbit_no;
    subcircuit.push_back(cnot);
    addCnfAdjacent(e.first.cbit_no,e.first.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }

  for(auto e:subcircuit){
    cout<<"subcircuit"<<endl;
    cout<<e.cbit_no<<endl;
  }

  //deleteSubcircuitgates(circuit,subcircuit);


}



void generateFirstSubcircuit(Graph &graph,vector<pair<Cnot,int>> &circuit,vector<vector <int>> &connect,vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,vector<Cnot> &subcircuit,int total_cnf_no,int cnf_fixed_no,int cnf_searched_answer_no,vector<string> &cnf_fixed,vector<string> &cnf_searched_answer,vector<string> &cnf_start,vector<int> vacancy_grid,state &start,vector<Cnot> &temp_subcircuit){

  int vertices_no_graph = circuit.size();
  vector<pair <Cnot,int> > temp_circuit;
  vector<int> result_int,progress;
  temp_circuit = circuit;
  for(auto e:circuit){
    subcircuit.push_back(e.first);
    addCnfAdjacent(e.first.cbit_no,e.first.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }
    //generateConstructedDependencyGraph(graph,circuit,connect,vertices_no_graph,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
    SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
    if(checkSAT()==1){
      cout<<"SWAP_NUMBER:0"<<endl;
    }
    while(subcircuit.size()!=0){
      deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
      subcircuit.pop_back();
    } 
    if(checkSAT()==-1){
      int success,fail;
      success=0;
      circuit = temp_circuit;
      fail = circuit.size();
      while(fail-success>1){
        circuit = temp_circuit;
        generateConstructedDependencyGraph(graph,circuit,connect,(success+fail)/2,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
        SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
        while(subcircuit.size()!=0){
          deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
          subcircuit.pop_back();
        }
        progress.push_back((success+fail)/2);
        if(checkSAT()==1) success = (success+fail)/2;
        else fail = (success+fail)/2;
      }
      circuit = temp_circuit;
      generateConstructedDependencyGraph(graph,circuit,connect,success,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
      SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
      cout<<"subcircuit"<<endl;
      for(auto e:subcircuit){
        cout<<e.cbit_no<<","<<e.tbit_no<<endl;
      }
      temp_subcircuit=subcircuit;
      while(subcircuit.size()!=0){
          deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
          subcircuit.pop_back();
      }
      deleteSubcircuitgates(circuit,temp_subcircuit);
      getCnfAnswer(result_int);
      generateState(start,node_list,result_int,vacancy_grid);
      cout<<"start placement"<<endl;
      printState(start);
      for(auto e:progress){
        cout<<e<<endl;
      }
      cout<<"success:"<<success<<endl;
    }

}

int generateFollowingSubcircuit(Graph &graph,vector<pair<Cnot,int>> &circuit,vector<vector <int>> &connect,vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,vector<Cnot> &subcircuit,int total_cnf_no,int cnf_fixed_no,int cnf_searched_answer_no,vector<string> &cnf_fixed,vector<string> &cnf_searched_answer,vector<string> &cnf_start,vector<int> vacancy_grid,state &start,int &success){

  int vertices_no_graph = circuit.size();
  vector< pair<Cnot,int> > temp_circuit;
  vector<Cnot> temp_subcircuit;
  vector<int> result_int,progress;
  int swap_no; 
  float swap_efficiency;
  vector< pair<vector<Cnot>,int> > answers,swap;
  state temp_start(w,h);
  temp_start = start;

  temp_circuit = circuit;

  for(auto e:circuit){
    subcircuit.push_back(e.first);
    addCnfAdjacent(e.first.cbit_no,e.first.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }
  total_cnf_no = cnf_fixed_no+cnf_adjacent_no+cnf_searched_answer_no;
  CnfInputGenerate(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no);
  SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
  if(checkSAT()==1) {
    pair<vector<Cnot>,int> answer;
    swap_no = Astar2(start,subcircuit);
    answer.second = swap_no;
    answer.first = subcircuit;
    answers.push_back(answer);
  }
  progress.push_back(circuit.size());
  while(subcircuit.size()!=0){
    deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
    subcircuit.pop_back();
  }
  cout<<"vertices"<<vertices_no_graph<<endl;
  cout<<subcircuit.size()<<endl;
  if(checkSAT()==1){
    getCnfAnswer(result_int);
    state temp_goal(w,h);
    generateState(temp_goal,node_list,result_int,vacancy_grid);
    swap_no = Astar(start,temp_goal);
    //swap_no = Astar2(start,subcircuit);
    cout<<"swap_no:"<<swap_no<<endl;
  }
  while(subcircuit.size()!=0){
     deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
     subcircuit.pop_back();
  } 
  int fail;
  success=0;
  circuit = temp_circuit;
  fail = circuit.size();
  int count = 0;
  while(fail-success>1){
    progress.push_back((success+fail)/2);
    circuit = temp_circuit;
    start = temp_start;
    generateConstructedDependencyGraph(graph,circuit,connect,(success+fail)/2,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
    total_cnf_no = cnf_fixed_no+cnf_adjacent_no+cnf_searched_answer_no;
    CnfInputGenerate(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no);
    SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
    if(checkSAT()==1){
      count++;
      success = (success+fail)/2;
      getCnfAnswer(result_int);
      pair<vector<Cnot>,int> answer;
      printState(start);
      for(auto e:subcircuit){
        cout<<e.cbit_no<<","<<e.tbit_no<<endl;
      }
      swap_no=Astar2(start,subcircuit);
      answer.first = subcircuit;
      answer.second = swap_no;
      answers.push_back(answer);
    }
    else fail = (success+fail)/2;
    while(subcircuit.size()!=0){
      deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
      subcircuit.pop_back();
    }
  }
  circuit = temp_circuit;
  start = temp_start;
  float max_efficiency=0;
  for(auto e:answers){
    if(e.second==0) swap_efficiency = 10000;
    else{
      swap_efficiency = float(e.first.size()/e.second);
      if(swap_efficiency>max_efficiency)
      max_efficiency = swap_efficiency;
    }
  }
  if(swap_efficiency == 10000){
    for(auto e:answers){
      int max_size = 0;
      if(e.second==0&&e.first.size()>max_size){
        temp_subcircuit=e.first;
        swap_no = e.second;
        max_size = temp_subcircuit.size();
      }
    }
  }
  else{
    for(auto e:answers){
      swap_efficiency = float(e.first.size()/e.second);
      cout<<"swap_efficiency"<<swap_efficiency<<endl;
      int max_size = 0;
      if(swap_efficiency==max_efficiency&&e.first.size()>max_size){
        temp_subcircuit=e.first;
        swap_no = e.second;
        max_size = temp_subcircuit.size();
      }
    }
  }

    start = temp_start;
    int a = Astar2(start,temp_subcircuit);
 
    deleteSubcircuitgates(circuit,temp_subcircuit);
    /*for(auto e:temp_subcircuit){
      addCnfAdjacent(e.cbit_no,e.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
    }
    total_cnf_no = cnf_fixed_no+cnf_adjacent_no+cnf_searched_answer_no;
    CnfInputGenerate(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no);
    SAT();*/
    cout<<"temp_subcircuit"<<endl;
    for(auto e:temp_subcircuit){
      cout<<e.cbit_no<<","<<e.tbit_no<<endl;
    }
    /* while(subcircuit.size()!=0){
        deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
        subcircuit.pop_back();
     }
    getCnfAnswer(result_int);
    generateState(start,node_list,result_int,vacancy_grid);*/
    
  cout<<"progress"<<endl;
  for(auto e:progress){
        cout<<e<<endl;
  }
  if(circuit.size()==temp_subcircuit.size()) circuit.clear();

  cout<<"swap"<<swap_no<<endl;

  return(success);

  
}

int generateFollowingSubcircuit2(Graph &graph,vector<pair<Cnot,int>> &circuit,vector<vector <int>> &connect,vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,vector<Cnot> &subcircuit,int total_cnf_no,int cnf_fixed_no,int cnf_searched_answer_no,vector<string> &cnf_fixed,vector<string> &cnf_searched_answer,vector<string> &cnf_start,vector<int> vacancy_grid,state &start,int &success,state &goal){

  int vertices_no_graph = circuit.size();
  vector< pair<Cnot,int> > temp_circuit;
  vector<int> progress,result_int;
  state temp_start(w,h);
  vector<Cnot> temp_subcircuit;
  temp_start = start;

  temp_circuit = circuit;
  temp_subcircuit.clear();
  subcircuit.clear();

  for(auto e:circuit){
    subcircuit.push_back(e.first);
    addCnfAdjacent(e.first.cbit_no,e.first.tbit_no,node_list,cnf_adjacent,cnf_adjacent_no);
  }
  SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
  if(checkSAT()==1){
     success = circuit.size();
     progress.push_back(circuit.size());
     temp_subcircuit = subcircuit;
     while(cnf_adjacent.size()!=0){
       deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
       subcircuit.pop_back();
    }
  }
  else{
    int fail = circuit.size();
    success=0;
    circuit = temp_circuit;
    while(fail-success>1){
      progress.push_back((success+fail)/2);
      circuit = temp_circuit;
      start = temp_start;
      generateConstructedDependencyGraph2(graph,circuit,connect,(success+fail)/2,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
      SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
      cout<<"subcircuit"<<endl;
      for(auto e:subcircuit){
        cout<<e.cbit_no<<","<<e.tbit_no<<endl;
      }
      for(auto e:cnf_adjacent){
        cout<<e<<endl;
      }
      if(checkSAT()==1) success = (success+fail)/2;
      else fail = (success+fail)/2;
      temp_subcircuit = subcircuit;
      while(cnf_adjacent.size()!=0){
        cout<<"delete"<<endl;
        deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
        subcircuit.pop_back();
        cout<<cnf_adjacent.size()<<endl;
       }
    }
  }

  generateConstructedDependencyGraph2(graph,circuit,connect,success,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
  SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
  
   
   circuit = temp_circuit;
   deleteSubcircuitgates(circuit,subcircuit);
   getCnfAnswer(result_int);
   generateState(goal,node_list,result_int,vacancy_grid);
   searchBestPlacement(result_int,cnf_start,cnf_fixed,cnf_adjacent,cnf_fixed_no,cnf_adjacent_no, vacancy_grid, start, goal, node_list);
   start = temp_start;

  //deleteSubcircuitgates(circuit,temp_subcircuit);
    
    cout<<"progress"<<endl;
    for(auto e:progress){
        cout<<e<<endl;
    }
    //if(circuit.size()==temp_subcircuit.size()) circuit.clear();

   return(success);

  
}

int searchSwap(Graph &graph,vector<pair<Cnot,int>> &circuit,vector<vector <int>> &connect,vector<Node> &node_list,vector<string> &cnf_adjacent,int &cnf_adjacent_no,vector<Cnot> &subcircuit,int total_cnf_no,int cnf_fixed_no,int cnf_searched_answer_no,vector<string> &cnf_fixed,vector<string> &cnf_searched_answer,vector<string> &cnf_start,vector<int> vacancy_grid,state &start,int &success,state &goal,vector<Cnot> &temp_fi_subcircuit){

  int vertices_no_graph = circuit.size();
  vector< pair<Cnot,int> > temp_circuit;
  vector<Cnot> temp_subcircuit;
  vector<int> result_int,progress,temp_result_int;
  int swap_no;
  float swap_efficiency;
  vector< pair<vector<Cnot>,int> > answers,swap;
  vector< pair<vector<int>, vector<Cnot> > > temp_states;
  state temp_start(w,h);
  temp_start = start;

  temp_circuit = circuit;

  circuit = temp_circuit;

  //generateDependencyGraph(graph,circuit,connect);


  int gate_number = success;
  /*if(success<8) gate_number = success;
  else gate_number = 15;*/

  while(gate_number<=success){
    int SAT_count = 0;
    state temp_goal(w,h);
    int temp_total_cnf_no=0;
    int state_eval;
    int temp_state_eval=100000;
    pair<vector<int>, vector<Cnot>> temp_state;
    vector<string> temp_cnf_start;
    vector<string> temp_cnf_adjacent;
    circuit = temp_circuit;
    start = temp_start;
    cnf_searched_answer.clear();
    cnf_searched_answer_no=0;
    generateConstructedDependencyGraph(graph,circuit,connect,gate_number,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit);
    cout<<"graph_subcircuit"<<endl;
    for(auto e:subcircuit){
      cout<<e.cbit_no<<","<<e.tbit_no<<endl;
    }
    SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
    pair<vector<Cnot>,int> answer;
    
    while(checkSAT()==1&&SAT_count<1000){

        getCnfAnswer(result_int);
        addCnfSearchedAnswer(result_int,cnf_searched_answer,cnf_searched_answer_no); 
        generateState(temp_goal,node_list,result_int,vacancy_grid);
        state_eval = manhattanDist(start,temp_goal) ;//+ differentNodeNo(start,temp_goal);
        //manhattan_dists.push_back(state_eval);
        if(state_eval<temp_state_eval){
          temp_state_eval = state_eval;
          goal = temp_goal;
          temp_result_int = result_int;
        }
        SAT(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no,cnf_fixed_no,cnf_adjacent_no,cnf_searched_answer_no);
        SAT_count++;
        std::cout<<SAT_count<<std::endl;
    }
    cout<<"start"<<endl;
    printState(start);
    cout<<"goal"<<endl;
    printState(goal);
    cout<<cnf_adjacent_no<<endl;
    swap_no = Astar(start,goal);
    answer.first = subcircuit;
    answer.second = swap_no;
    answers.push_back(answer);
    temp_state.first = temp_result_int;
    temp_state.second = subcircuit;
    temp_states.push_back(temp_state);

    while(cnf_adjacent.size()!=0){
      deleteCnfAdjacent(cnf_adjacent,cnf_adjacent_no);
      subcircuit.pop_back();
    }
    gate_number++;
  }
  circuit = temp_circuit;
  start = temp_start;
  float max_efficiency=0;
  for(auto e:answers){
    if(e.second==0) {
      swap_efficiency = 10000;
      break;
    }
    else{
      swap_efficiency = float(e.first.size()/e.second);
      if(swap_efficiency>max_efficiency)
      max_efficiency = swap_efficiency;
    }
  }


  if(swap_efficiency == 10000){
    for(auto e:answers){
      int max_size = 0;
      if(e.second==0&&e.first.size()>max_size){
        temp_subcircuit=e.first;
        swap_no = e.second;
        max_size = temp_subcircuit.size();
      }
      for(auto e:temp_states){
        if(isSameCircuit(e.second,temp_subcircuit)) result_int = e.first;  
    }
    }
  }
  else{
    for(auto e:answers){
      cout<<"adjacent_gate_no"<<e.first.size()<<" swap_no:"<<e.second<<endl;
      swap_efficiency = float(e.first.size()/e.second);
      cout<<e.first.size()<<"  "<<e.second<<endl;
      cout<<"swap_efficiency"<<float(e.first.size()/e.second)<<endl;
      int max_size = 0;
      if(swap_efficiency==max_efficiency&&e.first.size()>max_size){
        temp_subcircuit=e.first;
        swap_no = e.second;
        max_size = temp_subcircuit.size();
      }
    }
    for(auto e:temp_states){
        if(isSameCircuit(e.second,temp_subcircuit)) result_int = e.first;  
    }
  }

  for(auto e:answers){
    cout<<"la_adjacent_gate_no"<<e.first.size()<<" swap_no:"<<e.second<<endl;
  }

  generateState(goal,node_list,result_int,vacancy_grid);
  start = temp_start;
  swap_no = Astar(start,goal);
   
  temp_fi_subcircuit = temp_subcircuit; 
  deleteSubcircuitgates(circuit,temp_subcircuit);

  cout<<"swap"<<swap_no<<endl;

  return(swap_no);

}

auto main(int argc, char* argv[]) -> int{

   auto input = argc >= 2 ? argv[1] : "test.txt";

   int file_number = 1;
   int gate_number,experiment_number;

   cout<<"qubit number:";
   cin>>qubit_no;

   cout<<"gate number:";
   cin>>gate_number;
   cout<<"experiment number:";
   cin>>experiment_number;

  for(int i=0;i<experiment_number;i++){


   int suncircuit_no,grid_size,cnf_fixed_no,total_cnf_no,vacancy_qubit_no,total_swap_no,first_subcircuit,success;
   grid_size=cnf_fixed_no=vacancy_qubit_no=suncircuit_no=first_subcircuit=total_swap_no=total_cnf_no=success=0;
   vector<Node> node_list;
   vector<string> cnf_fixed;
   vector<int> subcircuits_no,manhattan_dists,vacancy_grid,swap_nos;

   first_subcircuit = 1;
   

   vector< pair<Cnot,int> > circuit;
   vector< vector <Cnot> > subcircuits;
   vector< pair<int,int> > adjacent_pairs;
   
   vector<Cnot> temp_final_subcircuit,subcircuit,temp_subcircuit;
   vector<int> result_swap_nos;

   inputFileGenerate(gate_number,qubit_no);

   circuitGenerate(circuit,"test.txt");
   /*for(auto e:circuit){
      gate_number++;
   }*/

   int swap_no,cnf_adjacent_no,node_dist,temp_node_dist,different_node_no,state_eval,cnf_searched_answer_no,temp_total_cnf_no;
   swap_no = cnf_adjacent_no = cnf_searched_answer_no = 0;
   vector<state> states;

   vector<string> cnf_start,cnf_adjacent,cnf_searched_answer;

   h=sqrt(qubit_no);
   if(h*h!=qubit_no) h++;
  
   w=qubit_no/h;
   if(h*w!=qubit_no) w++;

   grid_size=w*h;
   vacancy_qubit_no = grid_size - qubit_no;

   if(vacancy_qubit_no){
     int count = 1;
     for(int i=0;i<vacancy_qubit_no;i++){
       if(i%2==0){
         vacancy_grid.push_back(w*(h-1)+count);
       }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
       else {
         vacancy_grid.push_back(w*h-count+1);
         count++;
       }
      }
   }

   state start(w,h),goal(w,h),temp_goal(w,h);

   nodeGenerate(w,h,qubit_no,node_list);
   fixedCnfEncoder(node_list,cnf_fixed,cnf_fixed_no,vacancy_grid);

   Graph graph;
   vector<vector<int>> connect;
   connect.resize(2*gate_number);

   for(int i=0;i<2*gate_number;i++){
     connect[i].resize(2*gate_number);
   }

   generateFirstSubcircuit(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,temp_subcircuit);
   temp_goal = start;
   states.push_back(start);
   subcircuits.push_back(temp_subcircuit);
   printState(start);

   while(circuit.size()!=0){
      //for(int i=0;i<3;i++){
        cnf_adjacent.clear();
        success=generateFollowingSubcircuit2(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,success,goal);
        cout<<"success"<<success<<endl;
        //swap_no += Astar(start, goal);
        //result_swap_nos.push_back(swap_no);
        //swap_no = swap_no+searchSwap(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,success,goal,temp_fi_subcircuit);
        subcircuits.push_back(subcircuit);
        start = goal;
        states.push_back(goal);
        //swap_nos.push_back(swap_no);
   }

   //swap_no += Astar(start, temp_goal);

    states.push_back(temp_goal);

   cout<<circuit.size()<<endl;
   for(auto e:circuit){
     cout<<e.first.cbit_no<<","<<e.first.tbit_no<<","<<e.second<<endl;
   }

   for(auto e:states){
      printState(e);
      cout<<endl;
    }

    cout<<subcircuits.size()<<endl;
    for(auto e:subcircuits){
      for(auto f:e){
        cout<<f.cbit_no<<","<<f.tbit_no<<endl;
      }
      cout<<endl;
    }

    for(auto e:result_swap_nos){
      cout<<"swap_nos"<<endl;
      cout<<e<<endl;
    }

    /*for(auto e:evaluate_progress){
      cout<<"evaluates"<<endl;
      cout<<e<<endl;
    }*/

    cout<<"swap_no:"<<swap_no<<endl;

   benchmarkFileGenerate(states,file_number);
   file_number++;
    }



   return 0;
}


/*auto main(int argc, char* argv[]) -> int{

    auto input = argc >= 2 ? argv[1] : "test1.txt";

    int grid_size,gate_number,cnf_fixed_no,vacancy_qubit_no;
    grid_size=gate_number=cnf_fixed_no=vacancy_qubit_no=0;
    vector<Node> node_list;
    vector<string> cnf_fixed;
    vector<int> subcircuits_no,manhattan_dists,vacancy_grid,swap_nos;

    int suncircuit_no,first_subcircuit,total_swap_no,total_cnf_no,success,experiment_number;
    suncircuit_no=total_swap_no=total_cnf_no=0;
    first_subcircuit = 1;
    //std::vector<Cnot> circuit;
    std::vector<std::pair<Cnot,int>> circuit;
    vector< vector <Cnot> >subcircuits;
    std::vector<std::pair<int,int>> adjacent_pairs;
    std::vector<state> states;
    vector<Cnot> temp_fi_subcircuit,subcircuit,temp_subcircuit;
    vector<int> result_swap_nos;

    /*cout<<"qubit number:";
    cin>>qubit_no;

    cout<<"gate number:";
    cin>>gate_number;
     std::cout<<"experiment number:";
  std::cin>>experiment_number;

    qubit_no = circuitGenerate(circuit,input);
    for(auto e:circuit){
      gate_number++;
    }

    //for(int i=0;i<experiment_number;i++){

    

    int swap_no,cnf_adjacent_no,node_dist,temp_node_dist,different_node_no,state_eval,cnf_searched_answer_no,temp_total_cnf_no;
    swap_no = 0;

    std::vector<std::string> cnf_start,cnf_adjacent,cnf_searched_answer;
  
     h=std::sqrt(qubit_no);
     if(h*h!=qubit_no) h++;
  
     w=qubit_no/h;
     if(h*w!=qubit_no) w++;

     grid_size=w*h;
     vacancy_qubit_no = grid_size - qubit_no;
     
     if(vacancy_qubit_no){
       int count = 1;
       for(int i=0;i<vacancy_qubit_no;i++){
         if(i%2==0){
           vacancy_grid.push_back(w*(h-1)+count);
         }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
         else {
           vacancy_grid.push_back(w*h-count+1);
           count++;
         }
        }
     }

     inputFileGenerate(gate_number,grid_size);
     qubit_no=circuitGenerate(circuit,"test.txt");

     state start(w,h);
     state goal(w,h);
     state temp_goal(w,h);
     nodeGenerate(w,h,qubit_no,node_list);
     fixedCnfEncoder(node_list,cnf_fixed,cnf_fixed_no,vacancy_grid);

     Graph graph;
     vector<vector<int>> connect;
     connect.resize(2*gate_number);
     for(int i=0;i<2*gate_number;i++){
      connect[i].resize(2*gate_number);
    }

    generateFirstSubcircuit(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,temp_subcircuit);
    temp_goal = start;
    states.push_back(start);
    subcircuits.push_back(temp_subcircuit);
    while(circuit.size()!=0){
      //for(int i=0;i<5;i++){
        success=generateFollowingSubcircuit2(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,success);
        cout<<"success"<<success<<endl;
        swap_no = swap_no+searchSwap(graph,circuit,connect,node_list,cnf_adjacent,cnf_adjacent_no,subcircuit,total_cnf_no,cnf_fixed_no,cnf_searched_answer_no,cnf_fixed,cnf_searched_answer,cnf_start,vacancy_grid,start,success,goal,temp_fi_subcircuit);
        subcircuits.push_back(temp_fi_subcircuit);
        start = goal;
        states.push_back(goal);
        swap_nos.push_back(swap_no);
        
        //printState(start);   
    }

      states.push_back(temp_goal);

    int SAT_count = 0;
    int state_eval;
    int temp_state_eval=100000;
    vector<int> result_int,progress,temp_result_int;
    pair<vector<int>, vector<Cnot>> temp_state;
    vector<string> temp_cnf_start;
    vector<string> temp_cnf_adjacent;
    cnf_searched_answer.clear();
    cnf_searched_answer_no=0;
    total_cnf_no = cnf_fixed_no+cnf_adjacent_no+cnf_searched_answer_no;
    CnfInputGenerate(cnf_start,cnf_fixed,cnf_adjacent,cnf_searched_answer,total_cnf_no);
    SAT();

    while(checkSAT()==1&&SAT_count<200){

        getCnfAnswer(result_int);
        addCnfSearchedAnswer(result_int,cnf_searched_answer,cnf_searched_answer_no); 
        generateState(temp_goal,node_list,result_int,vacancy_grid);
        state_eval = manhattanDist(start,temp_goal) ;//+ differentNodeNo(start,temp_goal);
        //manhattan_dists.push_back(state_eval);
        if(state_eval<temp_state_eval){
          temp_state_eval = state_eval;
          goal = temp_goal;
          temp_result_int = result_int;
        }
        temp_total_cnf_no = total_cnf_no + cnf_searched_answer_no;
        CnfInputGenerate(temp_cnf_start,cnf_fixed,temp_cnf_adjacent,cnf_searched_answer,temp_total_cnf_no);
        SAT();
        SAT_count++;
        std::cout<<SAT_count<<std::endl;
    }
    printState(goal);
    swap_no += Astar(start,temp_goal);
    
    
    cout<<"circuit_size"<<circuit.size()<<endl;

    for(auto e:states){
      printState(e);
      cout<<endl;
    }

    cout<<subcircuits.size()<<endl;
    for(auto e:subcircuits){
      for(auto f:e){
        cout<<f.cbit_no<<","<<f.tbit_no<<endl;
      }
      cout<<endl;
    }

    for(auto e:swap_nos){
      cout<<e<<endl;
    }

    
    result_swap_nos.push_back(swap_no);
    //}

    total_swap_no=0;

    for(auto e:result_swap_nos){
      cout<<e<<endl;
      total_swap_no+=e;
    }

    

    cout<<"gates:"<<gate_number<<" "<<"qubit_no:"<<qubit_no<<endl;

    cout<<"result_swap_no"<<(total_swap_no/experiment_number)<<endl; 
    //std::ofstream file("test.dot");
    //boost::write_graphviz(file, graph);
     
    return 0;

}*/
