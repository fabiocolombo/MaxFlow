//Author: Fabio Colombo colombo.fabio@gmail.com

#ifndef MAXFLOW_HPP
#define MAXFLOW_HPP
#include <list>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <cassert>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>
#include "floatUtils.hpp"

#define ALPHA 6
#define BETA 12
#define FREQ 0.5;

namespace opt{

//Template class to solve Maximum Flow Problems
//using the push relabel algorithm.
//Based on:
//"On implementing the push—relabel method for the maximum flow problem"
//Lecture Notes in Computer Science, 1995, Volume 920/1995, 157-171
//by A. Goldberg and B. Cherkassky

template <class T> class MaxFlow{
private:
  /**node of the network*/
  struct fnode;
  /**arc of the network*/
  struct fedge{
    /**residual capacity*/
    T resCap;
    /**real capacity*/
    T realCap;
    /**pointer to the opposite arc
    fedge *rev;
    /**pointer to the tail node*/
    fnode *n1;
    /**pointer to the head node*/
    fnode *n2;
  };

  /**node of the network
  struct fnode{
    /**flow excess in the node*/
    T excess;
    /**node label*/
    int d;
    /**flag used by the breadth-first visit*/
    bool flag;
    /**node outcut*/
    std::list<fedge*> adj;
    /**node incut*/
    std::list<fedge*> inList;
    typename std::list<fedge*>::iterator current;
    /**pointer to the next node in the nodes list*/
    fnode *next;
    /**pointer to the previous node in the nodes list*/
    fnode *prev;
    /**index of this node in the array data*/
    int index;
  };

public:

  typedef std::map<boost::tuple<int,int>, boost::tuple<T,T,fedge*> > arcMap;

 private:
  /**nodes number*/
  int n;
  /**number of executed push moves*/
  int pushes;
  /**number of executed relabel moves*/
  int relabels;
  /**counter to handle the global update*/
  int relCount;
  /**source node pointer*/
  fnode *source;
  /**sink node pointer*/
  fnode *sink;
  /**source node index*/
  int so;
  /**sink node index*/
  int si;
  /**nodes vector*/
  std::vector<fnode*> nodes;
  /**data associated with the network arcs*/
  arcMap arcData;
  /**input arcs*/
  std::vector< boost::tuple<int,int> > realArc;
  /**head pointers to active nodes*/
  std::vector<fnode*> active;
  /**head pointers to inactive nodes*/
  std::vector<fnode*> inactive;

  /**biggest label between active nodes*/
  int aindex;
  /**biggest label between inactive nodes*/
  int iindex;

  /**create data structure needed by the algorithm*/
  void createDS();

  /**push move
     @param no node with excess flow
     @param e edge on which operate
  */
  void push(fedge *e, fnode *no);

  /**relabel move
     @param no node to be relabelled
  */
  void relabel(fnode *no);

  /**Discharge the node or deactivate it
     @param no node to be discharged
  */
  void discharge(fnode *no);

  /**Check if a node is active
     @param no node to check
     @return true if the node is active
  */
  bool isActive(fnode *no);
  /**Check if an arc is admissible
     @param from tail node
     @param to head node
  */
  bool admissible(fnode *from, fnode *to);

  /**Execute a global update*/
  void update();

  /**Add a node to an inactive nodes bucket
     @param no node to be added
     @param bucket inactive nodes bucket to use
  */
  void addInactive(fnode *no, int bucket);

  /**Remove a node from an inactive nodes bucket
     @param no node to be remove
     @param bucket inactive nodes bucket to use
  */
  void delInactive(fnode *no, int bucket);

  /**Add a node to an active nodes bucket
     @param no node to be added
     @param bucket active nodes bucket to use
  */
  void addActive(fnode *no, int bucket);

  /**Remove a node from an active nodes bucket
     @param no node to be remove
     @param bucket active nodes bucket to use
  */
  void delActive(fnode *no, int bucket);

  /**Reset data structure in order to restart the algorithm*/
  void reset();

  /**Greater than operator to specialized float precision version*/
  bool greater(const T& a, const T& b) const;

  /**return a const reference to the internal data*/
  const arcMap& getArcMap(){return arcData;}

public:

  /**Standard constructor
     @param n nodes number
  */
  MaxFlow(int n);

  /**Copy constructor
     @param mf solver to be copied
  */
  MaxFlow(const MaxFlow<T>& mf);

  ~MaxFlow();

  /**Add an arc to the network: must be invoked before solve() method
     @param from tail node
     @param to head node
     @param cap arc capacity
  */
  void addArc(int from, int to, T cap);

  /**Compute the minimum cut: must be invoked after solve() method
     @param soCut list in which will be inserted the arcs in the minimum cut near to the source
     @param siCut list in which will be inserted the arcs in the minimum cut near to the sink
  */
  void getCut(std::list< boost::tuple<int,int> >& soCut, std::list< boost::tuple<int,int> >& siCut);

  /**Change the capacity of an arc: to be invoked between two successive solve() invokation
     @param from tail node
     @param to new node
     @param cap new capacity
     @return true if has been possible to change the arc capacity
  */
  bool setCap(int from, int to, T cap);
  /**Select the source node: to be invoked before the solve() method
     @param source source node
  */
  void setSource(int source);
  /**Select the sink node: to be invoked before the solve() method
     @param sink sink node
  */
  void setSink(int sink);

  /**Compute the maximum flow value
     @return the maximum flow value
  */
  T solve();

  /**Print simple stats on the last algorithm execution*/
  void printStats();

  /**Compute the primal solution*/
  void secondStage();

  /**Compute the optimal flow on the an arc
     @param from tail node
     @param to head node
  */
  T getFlow(int from,int to);

  /**Check if the solution is feasible*/
  bool checkFlowAndCut();

  /**Print the problem in a compact form*/
  void printProblem(std::ostream& os=std::cout);
};

template <class T> MaxFlow<T>::MaxFlow(const MaxFlow<T>& mf):
n(mf.n), so(mf.so), si(mf.si),  nodes(n), arcData(mf.arcData){
  createDS();
}

template <class T> MaxFlow<T>::~MaxFlow(){
  for(typename std::vector<fnode*>::iterator it=nodes.begin();it!=nodes.end();++it)
    delete *it;
  if(!arcData.empty()){
	  typename arcMap::iterator it=arcData.begin();
          fedge *fPtr=boost::get<2>(it->second);
	  if(fPtr!=0){
		for(;it!=arcData.end();++it){
                  fPtr=boost::get<2>(it->second);
		  delete fPtr->rev;
		  delete fPtr;
		}
	  }
  }
}

template<class T> bool MaxFlow<T>::checkFlowAndCut(){
  for(int j=0;j<realArc.size();j++){
    if(getFlow(realArc[j].get<0>(),realArc[j].get<1>())<0) return false;
  }
  /**check flow constraints*/
  for(int i=0;i<n;++i){
    T af,fi=0;
    if(i==so) af=sink->excess;
    else if(i==si) af=-sink->excess;
    else af=0;
    for(int j=0;j<realArc.size();j++){
      if(realArc[j].get<0>()==i) fi+=getFlow(realArc[j].get<0>(),realArc[j].get<1>());
      else if(realArc[j].get<1>()==i) fi-=getFlow(realArc[j].get<0>(),realArc[j].get<1>());
    }
    if(af!=fi) return false;
  }
  return true;
}

template<class T> bool MaxFlow<T>::greater(const T& a, const T& b) const{
  return a>b;
}

template <class T> T MaxFlow<T>::getFlow(int from,int to){
  if(from>to){
    typename arcMap::iterator it=arcData.find(boost::make_tuple(to,from));
    if(it==arcData.end()) return -1;
    else return it->second.third->rev->realCap-it->second.third->rev->resCap;
  }
  else{
    typename arcMap::iterator it=arcData.find(boost::make_tuple(from,to));
    if(it==arcData.end()) return -1;
    else return it->second.third->realCap-it->second.third->resCap;
  }
}

template <class T> bool MaxFlow<T>::setCap(int from, int to, T cap){
  if(to<from){
    typename arcMap::iterator it=arcData.find(boost::make_tuple(to,from));
    if(it==arcData.end()) return false;
    boost::get<1>(it->second)=cap;
  }
  else{
    typename arcMap::iterator it=arcData.find(boost::make_tuple(from,to));
    if(it==arcData.end()) return false;
    boost::get<0>(it->second)=cap;
  }
  return true;
}

template <class T> void MaxFlow<T>::reset(){
  relabels=relCount=pushes=0;
  active.assign(active.size(),0);
  inactive.assign(inactive.size(),0);
  for(int i=0;i<n;++i){
    nodes[i]->next=0;
    nodes[i]->prev=0;
  }
  for(typename arcMap::iterator it=arcData.begin();it!=arcData.end();++it){
    fedge *e=boost::get<2>(it->second);
    e->resCap=boost::get<0>(it->second);
    e->rev->resCap=boost::get<1>(it->second);
  }
}


template <class T> MaxFlow<T>::MaxFlow(int n):n(n),pushes(0),relabels(0), so(0), si(0), nodes(n), active(n),inactive(n+1){}

template <class T> void MaxFlow<T>::setSource(int s){so=s;}

template <class T> void MaxFlow<T>::setSink(int s){si=s;}

template <class T> void MaxFlow<T>::secondStage(){
  const int BLACK=0;
  const int GREY=1;
  const int WHITE=2;
  //manage self-loops
  for(int i=0;i<n;++i){
    fnode *no= nodes[i];
    for(typename std::list<fedge*>::iterator ait=no->adj.begin();ait!=no->adj.end();++ait){
      if( (*ait)->n1 == (*ait)->n2 ){
        boost::tuple<T,T,fedge*> t=arcData[boost::make_tuple(i,i)];
        (*ait)->resCap=boost::get<1>(t);
      }
    }
  }

  //initialization
  for(int i=0;i<n;++i){
    nodes[i]->d=WHITE;
    active[i]=0;
    nodes[i]->current=nodes[i]->adj.begin();
  }
  fnode *ni,*nj,*restart,*tos,*bos,*r;
  fedge *a;
  T delta;
  tos=bos=0;

  //delete the cyclic flow and make a topologi nodes sort
  for(int i=0;i<n;++i){
    ni=nodes[i];
    if(ni->d==WHITE && greater(ni->excess,0) && ni!=source && ni!=sink){
      r=ni;
      r->d=GREY;
      do{
	for(;ni->current!=ni->adj.end();ni->current++){
	  a= *(ni->current);
	  if(a->realCap==0 && greater(a->resCap,0)){
	    nj=a->n2;
	    if(nj->d==WHITE){
	      nj->d=GREY;
	      addActive(ni,nj->index);
	      ni=nj;
	      break;
	    }
	    else
	      if(nj->d==GREY){
		/**find the minimum flow in the cycle*/
	        delta=a->resCap;
		while(true){
		  delta=std::min(delta,(*nj->current)->resCap);
		  if(nj==ni) break;
		  else nj=(*nj->current)->n2;
		}
		/**erase delta flow units from the flow*/
		nj=ni;
		while(true){
		  a=(*nj->current);
		  a->resCap-=delta;
		  a->rev->resCap+=delta;
		  nj=a->n2;
		  if(nj==ni) break;
		}
		/**backup DFS to the first saturated arc*/
		restart=ni;
		for(nj=(*ni->current)->n2; nj!=ni; nj=a->n2){
		  a= *(nj->current);
		  if(nj->d==WHITE || !greater(a->resCap,0)){
		    (*nj->current)->n2->d=WHITE;
		    if(nj->d!=WHITE) restart=nj;
		  }
		}
		if(restart!=ni){
		  ni=restart;
		  ++ni->current;
		  break;
		}
	      }
	  }
	}//end for

	if(ni->current==ni->adj.end()){
	  //finish the node i analysis
	  ni->d=BLACK;
	  if(ni!=source){
	    if(bos==NULL) bos=tos=ni;
	    else{
	      ni->next=tos;
	      tos=ni;
	    }
	  }
	  if(ni!=r){
	    ni=active[ni->index];
	    ni->current++;
	  }
	  else break;
	}
      }while(true);
    }

    /*return excesses*/
    /*note that sink is not on the stack*/
    if(bos!=0){
      typename std::list<fedge*>::iterator ait;
      for(ni=tos;ni!=bos;ni=ni->next){
	ait=ni->adj.begin();
	a= *ait;
	while(greater(ni->excess, 0)){
	  if(a->realCap==0 && greater(a->resCap,0)){
	    if(greater(ni->excess,a->resCap))
	      delta=a->resCap;
	    else
	      delta=ni->excess;
	    a->resCap-=delta;
	    a->rev->resCap+=delta;
	    ni->excess-=delta;
	    a->n2->excess+=delta;
	  }
	  ++ait;
	  if(ait==ni->adj.end()) break;
	  a=*ait;
	}
      }
      /*now do the bottom*/
      ni=bos;
      ait=ni->adj.begin();
      a=*ait;
      while(greater(ni->excess, 0)){
	if(a->realCap==0 && greater(a->resCap,0)){
	  if(greater(ni->excess,a->resCap))
	    delta=a->resCap;
	  else
	    delta=ni->excess;
	  a->resCap-=delta;
	  a->rev->resCap+=delta;
	  ni->excess-=delta;
	  a->n2->excess+=delta;
	}
	++ait;
	if(ait==ni->adj.end()) break;
	a=*ait;
      }
    }
  }
}

template <class T> void MaxFlow<T>::addArc(int from, int to, T cap){
  realArc.push_back(boost::make_tuple(from,to));
  if(to<from){
    typename arcMap::iterator it=arcData.find(boost::make_tuple(to,from));
	if(it==arcData.end()) arcData[boost::make_tuple(to,from)]=boost::make_tuple(T(0),cap,static_cast<fedge*>(0));
        else boost::get<1>(it->second)=cap;
  }
  else{
    typename arcMap::iterator it=arcData.find(boost::make_tuple(from,to));
    if(it==arcData.end()) arcData[boost::make_tuple(from,to)]=boost::make_tuple(cap,T(0),static_cast<fedge*>(0));
    else boost::get<0>(it->second)=cap;
  }
}

template <class T> void MaxFlow<T>::createDS(){
  for(int i=0; i<n;++i){
    nodes[i]=new fnode();
    nodes[i]->index=i;
  }
  source=nodes[so];
  sink=nodes[si];

  for(typename arcMap::iterator it=arcData.begin();it!=arcData.end();++it){
    fedge *e,*eRev;
    const boost::tuple<int,int>& p=it->first;
    boost::get<2>(it->second)=e=new fedge();
    eRev=new fedge();
    e->rev=eRev;
    eRev->rev=e;
    e->n1=nodes[p.get<0>()];
    e->n2=nodes[p.get<1>()];
    if(boost::get<0>(it->second)>0)
      e->n1->adj.push_front(e);
    else
      e->n1->adj.push_back(e);
    e->n2->inList.push_back(e);
    eRev->n1=nodes[p.get<1>()];
    eRev->n2=nodes[p.get<0>()];
    if(greater(boost::get<1>(it->second),0))
      eRev->n1->adj.push_front(eRev);
    else eRev->n1->adj.push_back(eRev);
    eRev->n2->inList.push_back(eRev);
    e->realCap=e->resCap=boost::get<0>(it->second);
    e->rev->realCap=e->rev->resCap=boost::get<1>(it->second);
  }
}



template <class T> T MaxFlow<T>::solve(){
  relCount=0;
  if(boost::get<2>(arcData.begin()->second)==0){
    createDS();
  }
  else reset();
  //reset node data
  for(int i=0;i<n;++i){
    nodes[i]->d=1;
    nodes[i]->excess=0;
    nodes[i]->flag=false;
  }
  source=nodes[so];
  sink=nodes[si];

  sink->d=0;
  source->d=n;


  //initialize preflow
  for(typename std::list<fedge*>::iterator it=source->adj.begin();it!=source->adj.end();++it){
    fedge *e=(*it);
    e->n2->excess+=e->resCap;
    e->rev->resCap=e->resCap;
    e->resCap=0;
    ++pushes;
  }

  //initialize buckets
  int lim=ALPHA*n+arcData.size();
  aindex=iindex=0;
  for(typename std::vector<fnode*>::iterator it=nodes.begin();it!=nodes.end();++it){
    (*it)->current=(*it)->adj.begin();
    if(*it!=source && *it!=sink ){
      int d=(*it)->d;
      if(isActive(*it)) addActive(*it,d);
      else addInactive(*it,d);
    }
  }

  while(aindex!=0){
    fnode *no=active[aindex];

    assert(aindex==no->d)      ;
    discharge(no);
    if( relCount*0.5 >lim){
      update();
      relCount=0;
    }
    else if(!greater(no->excess,0) && no->d<n){
      delActive(no,no->d);
      addInactive(no,no->d);
    }
  }
  return sink->excess;
}



template <class T> void MaxFlow<T>::push(fedge *e, fnode *no){
  ++pushes;
  T delta=std::min(no->excess, e->resCap);
  e->resCap-=delta;
  e->rev->resCap+=delta;
  e->n1->excess-=delta;
  e->n2->excess+=delta;
}

template <class T> void MaxFlow<T>::relabel(fnode *no){
  relCount+=BETA;
  ++relabels;
  int mind=n;
  for(typename std::list<fedge*>::iterator it=no->adj.begin();it!=no->adj.end();++it){
    relCount++;
    if(greater((*it)->resCap,0) && (*it)->n2->d < mind) mind=(*it)->n2->d;
  }
  no->d=mind+1;
}

template <class T> bool MaxFlow<T>::isActive(fnode *no){
  return greater(no->excess,0) && no->d<n;
}

template <class T> bool MaxFlow<T>::admissible(fnode *from, fnode *to){
  return from->d == to->d+1;
}

template <class T> void MaxFlow<T>::discharge(fnode *no){
  do{
    while(no->current != no->adj.end()){
      fedge *e = *(no->current);
      if(admissible(e->n1,e->n2) && greater(e->resCap,0) ){
	bool flag=(!greater(e->n2->excess,0));
	push(e,e->n1);
	if(e->n2!=source && e->n2!=sink && flag){
	  delInactive(e->n2,e->n2->d);
	  addActive(e->n2,e->n2->d);
	  e->n2->current=e->n2->adj.begin();
	}
	if(!greater(no->excess,0)) break;
      }
      no->current++;
    }
    if(greater(no->excess,0)){
      int oldd=no->d;
      delActive(no,no->d);
      relabel(no);
      if(active[oldd]==0 && inactive[oldd]==0){
	while(++oldd<n){
	  while(inactive[oldd]!=0){
	    inactive[oldd]->d=n;
	    delInactive(inactive[oldd],oldd);
	  }
	}
	break;
      }
      else if(no->d<n){
	addActive(no,no->d);
	no->current=no->adj.begin();
      }
      else break;
    }
    else{
      no->current=no->adj.begin();
      break;
    }
  }while(true);
}

template <class T> void MaxFlow<T>::getCut(std::list< boost::tuple<int,int> >& soCut, std::list< boost::tuple<int,int> >& siCut){
#ifdef PRINTDEBUG
  for(int i=0;i<n;++i)
    std::cout<<"nodo "<<i<<" livello:"<<nodes[i]->d<<" eccesso:"<<nodes[i]->excess<<std::endl;
#endif

  for(std::vector< boost::tuple<int,int> >::iterator it=realArc.begin();it!=realArc.end();++it){
    if( nodes[it->get<0>()]->d>=n && nodes[it->get<1>()]->d<n  ){
      soCut.push_back(*it);
    }
  }

  //get tside
  for(int i=0;i<n;++i) nodes[i]->flag=false;
  sink->flag=true;
  std::list<fnode*> queue;
  std::set<fnode*> tside;
  tside.insert(sink);
  queue.push_back(sink);
  do{
    fnode *v=queue.front();
    queue.pop_front();
    for(typename std::list<fedge*>::iterator it=v->inList.begin();it!=v->inList.end();++it){
      fedge *e=(*it);
      if(greater(e->resCap,0) && e->n1->flag==false && e->n1!=source){
	e->n1->flag=true;
	tside.insert(e->n1);
	queue.push_back(e->n1);
      }
    }
  }while(!queue.empty());
  //get siCut
  for(std::vector< boost::tuple<int,int> >::iterator it=realArc.begin();it!=realArc.end();++it){
    if(tside.find(nodes[it->get<0>()])==tside.end() && tside.find(nodes[it->get<1>()])!=tside.end())
      siCut.push_back(*it);
  }

}

template <class T> void MaxFlow<T>::update(){
  for(int i=0;i<n;++i){
    nodes[i]->flag=false;
    nodes[i]->d=n;
    active[i]=0;
    inactive[i]=0;
  }
  aindex=iindex=0;
  sink->flag=true;
  sink->d=0;
  std::list<fnode*> queue;
  queue.push_back(sink);
  do{
    fnode *v=queue.front();
    queue.pop_front();
    for(typename std::list<fedge*>::iterator it=v->inList.begin();it!=v->inList.end();++it){
      fedge *e= (*it);
      if(greater(e->resCap,0) && e->n1->flag==false && e->n1!=source){
	e->n1->d=std::min(n,e->n2->d+1);
	if(greater(e->n1->excess,0)){
	  addActive(e->n1,e->n1->d);
	}
	else{
	  addInactive(e->n1,e->n1->d);
	}
      	e->n1->flag=true;
	queue.push_back(e->n1);
      }
    }
  }while(!queue.empty());
  source->d=n;
}

template <class T> void MaxFlow<T>::printStats(){
  std::cout<<"push:"<<pushes<<std::endl;
  std::cout<<"relabel:"<<relabels<<std::endl;
}

template <class T> void MaxFlow<T>::addActive(fnode *no, int bucket){
  if(active[bucket]==0){
    no->next=no->prev=0;
    active[bucket]=no;
    if(bucket>aindex) aindex=bucket;
  }
  else{
    no->next=active[bucket];
    no->next->prev=no;
    no->prev=0;
    active[bucket]=no;
  }
}

template <class T> void MaxFlow<T>::delActive(fnode *no, int bucket){
  if(no->prev==0){
    active[bucket]=no->next;
    if(no->next!=0) no->next->prev=0;
    else{
      while(--aindex>0 && active[aindex]==0);
    }
  }
  else{
    no->prev->next=no->next;
    if(no->next!=0) no->next->prev=no->prev;
  }
  no->next=no->prev=0;
}

template <class T> void MaxFlow<T>::addInactive(fnode *no, int bucket){
  if(inactive[bucket]==0){
    no->next=no->prev=0;
    inactive[bucket]=no;
    if(bucket>iindex) iindex=bucket;
  }
  else{
    no->next=inactive[bucket];
    no->next->prev=no;
    no->prev=0;
    inactive[bucket]=no;
  }
}

template <class T> void MaxFlow<T>::delInactive(fnode *no,int bucket){
  if(no->prev==0){
    inactive[bucket]=no->next;
    if(no->next!=0) no->next->prev=0;
    else{
      while(--iindex>0 && inactive[iindex]==0);
    }
  }
  else{
    no->prev->next=no->next;
    if(no->next!=0) no->next->prev=no->prev;
  }
  no->next=no->prev=0;
}

template <class T> void MaxFlow<T>::printProblem(std::ostream& os){
  os<<"p max "<<n<<" "<<realArc.size()<<std::endl;
  os<<"n "<<so<<" s"<<std::endl;
  os<<"n "<<si<<" t"<<std::endl;
  for(std::vector< boost::tuple<int,int> >::iterator it=realArc.begin();it!=realArc.end();++it){
    if(it->get<1>()<it->get<0>()){
      typename arcMap::iterator ait=arcData.find(boost::make_tuple(it->get<1>(),it->get<0>()));

      os<<"a "<<it->get<0>()<<" "<<it->get<1>()<<" "<<boost::get<1>(ait->second)<<std::endl;
    }
    else{
      typename arcMap::iterator ait=arcData.find(boost::make_tuple(it->get<0>(),it->get<1>()));
      os<<"a "<<it->get<0>()<<" "<<it->get<1>()<<" "<<boost::get<0>(ait->second)<<std::endl;
    }
  }
}

template <> bool MaxFlow<double>::greater(const double& a, const double& b) const{
	return util::floatUtils::gr(a,b,1e-4);
}

}//end namespace

#endif //MAXFLOW_HPP
