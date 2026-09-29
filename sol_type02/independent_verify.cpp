// Independent enumeration: scan all 2^21 tournaments, then find group orbits
// by breadth-first closure under generators (not Python's backtracking or
// full-permutation orbit generation).  Standard C++17, no external solver.
#include <array>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <iostream>
#include <chrono>
#include <cassert>
#include <cstdint>
using namespace std;
using Perm=array<int,7>;
vector<pair<int,int>> edges;
int idx[7][7];
uint32_t act(uint32_t m,const Perm&p,bool reverse){
  bool a[7][7]{};
  for(unsigned e=0;e<edges.size();++e){
    auto [i,j]=edges[e];
    if((m>>e)&1u) a[p[i]][p[j]]=true;
    else a[p[j]][p[i]]=true;
  }
  uint32_t q=0;
  for(unsigned e=0;e<edges.size();++e){auto [i,j]=edges[e];
    bool bit=a[i][j];if(reverse)bit=!bit;if(bit)q|=1u<<e;
  }
  return q;
}
Perm identity(){return {0,1,2,3,4,5,6};}
pair<Perm,bool> sw(int a,int b){Perm p=identity();swap(p[a],p[b]);return {p,false};}
int main(){
 auto start=chrono::steady_clock::now();
 for(int i=0;i<7;i++)for(int j=i+1;j<7;j++){idx[i][j]=edges.size();edges.emplace_back(i,j);}
 array<array<int,7>,3> target={array<int,7>{3,3,3,3,3,3,3},
   array<int,7>{3,2,4,3,3,3,3},array<int,7>{2,2,4,4,3,3,3}};
 array<set<uint32_t>,3> valid;
 for(uint32_t mask=0;mask<(1u<<21);++mask){
   array<int,7>d{};
   for(unsigned e=0;e<21;++e){auto[i,j]=edges[e];++d[((mask>>e)&1)?i:j];}
   for(int t=0;t<3;++t)if(d==target[t])valid[t].insert(mask);
 }
 vector<pair<Perm,bool>> gen[3];
 for(int i=0;i<6;i++)gen[0].push_back(sw(i,i+1));
 gen[0].push_back({identity(),true});
 for(int i=3;i<6;i++)gen[1].push_back(sw(i,i+1));
 auto p=identity();swap(p[1],p[2]);gen[1].push_back({p,true});
 gen[2]={sw(0,1),sw(2,3),sw(4,5),sw(5,6)};
 p=identity();swap(p[0],p[2]);swap(p[1],p[3]);gen[2].push_back({p,true});
 const char*names[3]={"type0","intersecting","disjoint"};
 cout<<"{\n  \"boards_scanned\": "<<(1u<<21)<<",\n";
 for(int t=0;t<3;t++){
   set<uint32_t>unseen=valid[t];map<size_t,size_t>hist;vector<uint32_t>reps;
   while(!unseen.empty()){
     uint32_t seed=*unseen.begin();reps.push_back(seed);set<uint32_t>orb{seed};queue<uint32_t>q;q.push(seed);
     while(!q.empty()){
       auto x=q.front();q.pop();
       for(auto [p,r]:gen[t]){
         auto y=act(x,p,r);assert(valid[t].count(y));
         if(orb.insert(y).second)q.push(y);
       }
     }
     for(auto x:orb){assert(unseen.erase(x)==1);}hist[orb.size()]++;
   }
   cout<<"  \""<<names[t]<<"\": {\"labelled_count\": "<<valid[t].size()<<", \"orbit_count\": "<<reps.size()<<", \"representative_masks\": [";
   for(size_t k=0;k<reps.size();k++){if(k)cout<<", ";cout<<reps[k];}cout<<"], \"orbit_size_histogram\": {";
   bool first=true;for(auto [sz,n]:hist){if(!first)cout<<", ";first=false;cout<<"\""<<sz<<"\": "<<n;}cout<<"}},\n";
 }
 double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
 cout<<"  \"elapsed_seconds\": "<<seconds<<"\n}\n";
}
