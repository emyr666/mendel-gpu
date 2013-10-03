#include<iostream>
using namespace std;

class parent{
public:
  parent(int i);
  virtual ~parent();
  void init();
  virtual void init3(int)=0;
protected:
  virtual void init2()=0;
};

class child:public parent{
public:
  child(int i);
  ~child();
  void init3(int);
private:
  void init2();
};

parent::parent(int i){
  cout<<"Parent constructor with child val "<<i<<"\n";
}

parent::~parent(){
  cout<<"Parent destructor\n";
}

void parent::init3(int a){
  cout<<"Parent init3\n";
}

void child::init3(int a){
  parent::init3(a);
  cout<<"Child init3\n";
}

child::child(int i):parent(i){
  cout<<"Child constructor\n";
}

child::~child(){
  cout<<"Child destructor\n";
}

void parent::init(){
  cout<<"Parent init\n";
  init2();
}

void child::init2(){
  cout<<"Subclass init\n";
}

int main(){
  parent * p = new child(5);
  p->init();
  p->init3(5);
  delete p;


}
