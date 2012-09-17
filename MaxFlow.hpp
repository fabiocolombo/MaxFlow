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


/**classe per la gestione dei problemi
   di flusso massimo e taglio minimo,
   basata sull'implementazione di Andrew Goldberg
   dell'algoritmo push_relabel con global update
   e gap heuristics
*/
template <class T> class MaxFlow{
private:
  struct fnode;
  /**arco della rete di flusso*/
  struct fedge{
    /**capacit� residua*/
    T resCap;
    /**capacit� reale*/
    T realCap;
    /**puntatore all'arco opposto*/
    fedge *rev;
    /**puntatore al nodo di partenza*/
    fnode *n1;
    /**puntatore al nodo di arrivo*/
    fnode *n2;
  };
  /**nodo della rete di flusso*/
  struct fnode{
    /**eccesso associato al nodo*/
    T excess;
    /**label associata al nodo*/
    int d;
    /**flag utilizzata per visita in ampiezza*/
    bool flag;
    /**lista degli archi uscenti*/
    std::list<fedge*> adj;
    /**lista degli archi entranti*/
    std::list<fedge*> inList;
    typename std::list<fedge*>::iterator current;
    /**successore nella lista dei nodi attivi o in
       quella dei nodi inattivi*/
    fnode *next;
    /**predecessore nella lista dei nodi attivi o in
       quella dei nodi inattivi*/
    fnode *prev;
    /**indice nel vettore*/
    int index;
  };
public:
  typedef std::map<boost::tuple<int,int>, boost::tuple<T,T,fedge*> > arcMap;
 private:
  /**numero di nodi*/
  int n;
  /**numero di push effettuati*/
  int pushes;
  /**numero di relabel effettuati*/
  int relabels;
  /**numero che determina se avviare la global update*/
  int relCount;
  /**nodo sorgente*/
  fnode *source;
  /**nodo pozzo*/
  fnode *sink;
  /**indice sorgente*/
  int so;
  /**indice pozzo*/
  int si;
  /**vettore dei nodi della rete di flusso*/
  std::vector<fnode*> nodes;
  /**dati associati agli archi della rete di flusso*/
  arcMap arcData;
  /**archi del problema*/
  std::vector< boost::tuple<int,int> > realArc;
  /**teste delle liste dei nodi attivi*/
  std::vector<fnode*> active;
  /**teste delle liste dei nodi inattivi*/
  std::vector<fnode*> inactive;

  /**label pi� grande di un nodo attivo*/
  int aindex;
  /**label pi� grande di un nodo inattivo*/
  int iindex;

  /**crea le strutture dati utilizzate dall'algoritmo*/
  void createDS();
  /**operazione push
     @param no nodo con flusso in eccesso
     @param e lato su cui effettuare l'operazione
  */
  void push(fedge *e, fnode *no);
  /**operazione relabel
     @param no nodo da rietichettare
  */
  void relabel(fnode *no);
  /**Svuta un nodo o lo disattiva
     @param no nodo da svuotare
  */
  void discharge(fnode *no);
  /**Verifica se un nodo � attivo
     @param no nodo da verificare
     @return true se il nodo � attivo
  */
  bool isActive(fnode *no);
  /**Verifica se un certo arco � ammissibile
     @param from nodo di partenza
     @param to nodo di arrivo
  */
  bool admissible(fnode *from, fnode *to);
  /**Global update*/
  void update();
  /**Aggiunge un nodo alla lista dei nodi inattivi di un certo livello
     @param no nodo da aggiungere
     @param bucket livello a cui aggiungere il nodo
  */
  void addInactive(fnode *no, int bucket);
  /**Cancella un nodo alla lista dei nodi inattivi di un certo livello
     @param no nodo da aggiungere
     @param bucket livello a cui aggiungere il nodo
  */
  void delInactive(fnode *no, int bucket);
  /**Aggiunge un nodo alla lista dei nodi attivi di un certo livello
     @param no nodo da aggiungere
     @param bucket livello a cui aggiungere il nodo
  */
  void addActive(fnode *no, int bucket);
  /**Cancella un nodo alla lista dei nodi attivi di un certo livello
     @param no nodo da aggiungere
     @param bucket livello a cui aggiungere il nodo
  */
  void delActive(fnode *no, int bucket);
  /**Resetta le strutture dati in modo da riavviare poi l'algoritmo*/
  void reset();
  /**Operazione di > utilizzata per specializzazione template double*/
  bool greater(const T& a, const T& b) const;
public:
  /**Costruisce un MaxFlow<T>
     @param n numero nodi del problema
  */
  MaxFlow(int n);
  /**Costruttore per copia
     @param mf oggetto da copiare
  */
  MaxFlow(const MaxFlow<T>& mf);
  /**Distruttore*/
  ~MaxFlow();
  /**Aggiunge un arco del problema, da invocare prima
     della chiamata a solve
     @param from nodo di partenza
     @param to nodo di arrivo
     @param cap capacit� dell'arco
  */
  void addArc(int from, int to, T cap);
  /**Calcola il taglio minimo: da chiamare dopo la chiamata a solve
     @param soCut lista in cui verranno inseriti i lati del taglio pi� vicino alla sorgente
     @param siCut lista in cui verranno inseriti i lati del taglio pi� vicino al pozzo
  */
  void getCut(std::list< boost::tuple<int,int> >& soCut, std::list< boost::tuple<int,int> >& siCut);
  /**Modifica la capacit� di un arco, da chiamare tra due
     chiamate di solve successive
     @param from nodo di partenza dell'arco
     @param to nodo di arrivo dell'arco
     @param cap nuova capacita
     @return true se � stato possibile modificare la capacit�
  */
  bool setCap(int from, int to, T cap);
  /**Seleziona/Modifica il nodo sorgente da chiamare prima di una chiamata
     a solve
     @param source nodo sorgente
  */
  void setSource(int source);
  /**Seleziona/Modifica il nodo destinazione da chiamare prima di una chiamata
     a solve
     @param sink nodo destinazione
  */
  void setSink(int sink);
  /**Calcola il valore del flusso massimo
     @return valore del flusso massimo
  */
  T solve();
  /**Stampa su cout statistiche dell'esecuzione dell'algoritmo*/
  void printStats();
  /**Costruisce la soluzione primale*/
  void secondStage();
  /**Ottiene flusso lungo un arco*/
  T getFlow(int from,int to);
  /**controlla se la soluzione � un flusso ammissibile*/
  bool checkFlowAndCut();
  /**ritorna un riferimento costante alla mappa interna degli archi*/
  const arcMap& getArcMap(){return arcData;}

  int getIndex(std::vector<fnode*>& v, fnode* no){
    return no->index;
  }

  /**stampa il problema*/
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
  /**verifico vincoli flusso*/
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
  //gestione cappi
  for(int i=0;i<n;++i){
    fnode *no= nodes[i];
    for(typename std::list<fedge*>::iterator ait=no->adj.begin();ait!=no->adj.end();++ait){
      if( (*ait)->n1 == (*ait)->n2 ){
        boost::tuple<T,T,fedge*> t=arcData[boost::make_tuple(i,i)];
        (*ait)->resCap=boost::get<1>(t);
      }
    }
  }
  //inizializzazione
  for(int i=0;i<n;++i){
    nodes[i]->d=WHITE;
    active[i]=0;
    nodes[i]->current=nodes[i]->adj.begin();
  }
  fnode *ni,*nj,*restart,*tos,*bos,*r;
  fedge *a;
  T delta;
  tos=bos=0;
  //elimina il flusso ciclico e ordina topologicamente i nodi
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
		/**cerca il flusso minimo nel ciclo*/
	        delta=a->resCap;
		while(true){
		  delta=std::min(delta,(*nj->current)->resCap);
		  if(nj==ni) break;
		  else nj=(*nj->current)->n2;
		}
		/**rimuove delta unita di flusso dal ciclo*/
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
	  //ho finito la scansione di i
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
    //std::cout<<"next--->"<<nodes[i]->next<<" prev--->"<<nodes[i]->prev<<std::endl;
    //std::cout<<"active--->"<<active[i]<<" inactive--->"<<inactive[i]<<std::endl;
  }
  source=nodes[so];
  sink=nodes[si];

  sink->d=0;
  source->d=n;

//  std::cout<<"compute flow from "<<so<<" to "<<si<<std::endl;
//  typedef std::map<boost::tuple<int,int>, boost::tuple<T,T,fedge*> > arcMap;
//
//  std::cout<<"arcData.size()="<<arcData.size()<<std::endl;
//  for(typename arcMap::iterator it=arcData.begin();it!=arcData.end();++it){
//      boost::tuple<T,T,fedge*> t=it->second;
//      fedge& e= *boost::get<2>(t);
//      if(e.realCap > 0)
//        std::cout<<"("<<e.n1->index<<","<<e.n2->index<<"):"<<e.realCap<<std::endl;
//      if(e.rev->realCap > 0)
//        std::cout<<"("<<e.rev->n1->index<<","<<e.rev->n2->index<<"):"<<e.rev->realCap<<std::endl;
//  }


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
    //    std::cout<<"aindex--->"<<aindex<<std::endl;
    #ifdef PRINTDEBUG
    std::cout<<"_________________update"<<std::endl;
    for(int i=0;i<n;++i){
     std::cout<<"nodo "<<i<<" "<<" :"<<nodes[i]->d<<": "<<nodes[i]->excess;
     if(isActive(nodes[i])) std::cout<<" A ";
     std::cout<<std::endl;
    }
    std::cout<<"archi:"<<std::endl;
    for(typename arcMap::iterator it=arcData.begin();it!=arcData.end();++it){
      fedge *e=it->second.get<2>();
      std::cout<<"analizzo il lato ("<<getIndex(nodes,e->n1)<<","<<getIndex(nodes,e->n2)<<") "<<"resCap="<<e->resCap<<std::endl;
      e=e->rev;
      std::cout<<"analizzo il lato ("<<getIndex(nodes,e->n1)<<","<<getIndex(nodes,e->n2)<<") "<<"resCap="<<e->resCap<<std::endl;
    }
    #endif
    fnode *no=active[aindex];
    //std::cout<<"scarico nodo "<<no<<std::endl;
    #ifdef PRINTDEBUG
    std::cout<<"lavoro su nodo: "<<getIndex(nodes,no)<<" eccesso:"<<no->excess<<std::endl;
    #endif
    // if(true){
    //   fnode *ptr=active[2];
    //   std::cout<<"active[aindex]:";
    //   if(ptr==0) std::cout<<"[]"<<std::endl;
    //   else{
    // 	std::cout<<"["<<ptr->index<<":"<<ptr->d;
    // 	ptr=ptr->next;
    // 	while(ptr!=0){
    // 	  std::cout<<", "<<ptr->index<<":"<<ptr->d;
    // 	  ptr=ptr->next;
    // 	}
    // 	std::cout<<"]"<<std::endl;
    //   }
    //   std::cout<<"active[no->d]:";
    //   ptr=active[no->d];
    //   if(ptr==0) std::cout<<"[]"<<std::endl;
    //   else{
    // 	std::cout<<"["<<ptr->index<<":"<<ptr->d;
    // 	ptr=ptr->next;
    // 	while(ptr!=0){
    // 	  std::cout<<", "<<ptr->index<<":"<<ptr->d;
    // 	  ptr=ptr->next;
    // 	}
    // 	std::cout<<"]"<<std::endl;
    //   }
    //   std::cout<<"aindex="<<aindex<<std::endl;
    //   std::cout<<"no->d="<<no->d<<std::endl;
    //   std::cout<<"no->index="<<no->index<<std::endl;
    // }
    assert(aindex==no->d)      ;
    discharge(no);
    #ifdef PRINTDEBUG
    //std::cout<<"eccesso dopo scarico:"<<no->excess<<std::endl;
    #endif
    if( relCount*0.5 >lim){
      #ifdef PRINTDEBUG
      std::cout<<"global update"<<std::endl;
      #endif
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
	//std::cout<<"eccesso prima:"<<e->n2->excess<<std::endl;
	push(e,e->n1);
	if(e->n2!=source && e->n2!=sink && flag){
	  //std::cout<<"il nodo "<<e->n2<<" diventa attivo"<<std::endl;
	  //std::cout<<"eccesso:"<<e->n2->excess<<std::endl;
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
      //      std::cout<<"faccio relabel di "<<no<<std::endl;
      delActive(no,no->d);
      relabel(no);
      //      std::cout<<"nuovo label:"<<no->d<<std::endl;
      if(active[oldd]==0 && inactive[oldd]==0){
	/* if(no->index==88){ */
	/*   std::cout<<"Lista bucket attivo "<<2<<":"<<std::endl; */
	/*   fnode *ptr=active[2]; */
	/*   if(ptr==0) std::cout<<"[]"<<std::endl; */
	/*   else{ */
	/*     std::cout<<"["<<ptr->index<<":"<<ptr->d; */
	/*     ptr=ptr->next; */
	/*     while(ptr!=0){ */
	/*       std::cout<<", "<<ptr->index<<":"<<ptr->d; */
	/*       ptr=ptr->next; */
	/*     } */
	/*     std::cout<<"]"<<std::endl; */
	/*   } */
	/*   std::cout<<"vecchio livello nodo 88--->"<<no->d<<std::endl; */
	/* } */
	no->d=n;
	/* if(no->index==88){ */
	/*   if(active[2]!=0) std::cout<<"active[2]->index:"<<active[2]->index<<std::endl; */
	/*   std::cout<<"setto 88 a"<<n<<std::endl; */
	/* } */
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

//#define DEB

template <class T> void MaxFlow<T>::addActive(fnode *no, int bucket){
#ifdef DEB
  if(no->index==88){
    std::cout<<"Lista bucket attivo "<<bucket<<":"<<std::endl;
    fnode *ptr=active[bucket];
    if(ptr==0) std::cout<<"[]"<<std::endl;
    else{
      std::cout<<"["<<ptr->index<<":"<<ptr->d;
      ptr=ptr->next;
      while(ptr!=0){
	std::cout<<", "<<ptr->index<<":"<<ptr->d;
	ptr=ptr->next;
      }
      std::cout<<"]"<<std::endl;
    }
  // assert(bucket==no->d);
  // std::cout<<"aggiungo il nodo "<<no<<" al bucket "<<bucket<<std::endl;
  // std::cout<<"no->next:"<<no->next<<" no->prev:"<<no->prev<<std::endl;
  // std::cout<<"testa lista:"<<active[bucket]<<std::endl;
    std::cout<<"aggiungo al bucket attivo "<<bucket<<", il nodo "<<no->index<<std::endl;
  }
#endif
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
#ifdef DEB
  if(no->index==88){
    std::cout<<"Lista bucket attivo dopo "<<bucket<<":"<<std::endl;
    fnode *ptr=active[bucket];
    if(ptr==0) std::cout<<"[]"<<std::endl;
    else{
      std::cout<<"["<<ptr->index<<":"<<ptr->d;
      ptr=ptr->next;
      while(ptr!=0){
	std::cout<<", "<<ptr->index<<":"<<ptr->d;
	ptr=ptr->next;
      }
      std::cout<<"]"<<std::endl;
    }
  }
#endif
}

template <class T> void MaxFlow<T>::delActive(fnode *no, int bucket){
#ifdef DEB
  if(no->index==88){
    std::cout<<"Lista bucket attivo "<<bucket<<":"<<std::endl;
    fnode *ptr=active[bucket];
    if(ptr==0) std::cout<<"[]"<<std::endl;
    else{
      std::cout<<"["<<ptr->index<<":"<<ptr->d;
      ptr=ptr->next;
      while(ptr!=0){
	std::cout<<", "<<ptr->index<<":"<<ptr->d;
	ptr=ptr->next;
      }
      std::cout<<"]"<<std::endl;
    }
    std::cout<<"cancello il nodo "<<no<<" dal bucket "<<bucket<<std::endl;
  }
#endif
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
#ifdef DEB
  if(no->index==88){
    std::cout<<"Lista bucket attivo dopo "<<bucket<<":"<<std::endl;
    fnode *ptr=active[bucket];
    if(ptr==0) std::cout<<"[]"<<std::endl;
    else{
      std::cout<<"["<<ptr->index<<":"<<ptr->d;
      ptr=ptr->next;
      while(ptr!=0){
	std::cout<<", "<<ptr->index<<":"<<ptr->d;
	ptr=ptr->next;
      }
      std::cout<<"]"<<std::endl;
    }
  }
#endif
  // if(bucket==1 && active[1]!=0){
  //   fnode *ptr=active[1];
  //   assert(active[1]->next!=active[1]);
  // }
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
