#include "slist.h"
#include <iostream>
#include <string>

#define quote(x) #x

using namespace std;

/*

Class Library File

*/


// Constructor
Slist::Slist() {
	sizeInt = 0;
	head = new CustomNode;
	tail = new CustomNode;
	tail->ptr = head;
}

// Destructor
Slist::~Slist() {
	cout << "Destructor" << endl;
	delete tail;
	clear();
	delete head;
}

// add(value)				//Adds a new value to the end of this list.
void Slist::add(CustomNode* value) {
	value->ptr = nullptr;
	//cout << "Tee" << endl;
	if(sizeInt == 0) {
		head->ptr = value;
	} else {
		tail->ptr->ptr = value;
	}
	tail->ptr = value;
	sizeInt++;
}

// clear()					//Removes all elements from this list.
void Slist::clear() {
	for(int i = sizeInt; i>=0; --i) {
		delete get(i);
	}
	head = nullptr;
	sizeInt = 0;
}

// equals(list)				//Returns true if the two lists contain the same elements in the same order.
bool Slist::equals(Slist* list) {
	if(list->size() != sizeInt) return false;
	for(int i = 0; i<sizeInt; ++i) {
		//cout << get(i)->printData() << " " << list->get(i)->toString() << endl;
		if(get(i)->toString() != list->get(i)->toString()) return false;
	}
	return true;
}

//get(index)				//Returns the element at the specified index in this list.
CustomNode* Slist::get(int index) {
	if(index == 0) {
		return head->ptr;
	}
	int ind = 0;
	CustomNode* next = head;
	while(ind < index && next->ptr != nullptr) {
		next = next->ptr;
		ind++;
	}
	return next->ptr;
}

//insert(index, value)		//Inserts the element into this list before the specified index.
void Slist::insert(int index, CustomNode* value) {
	CustomNode* before;
	if(index == 0) before = head;
	else before = get(index-1);
	CustomNode* after = before->ptr;
	before->ptr = value;
	value->ptr = after;
	sizeInt++;
}

//exchg(index1, index2)		//Switches the payload data of specified indexes
void Slist::exchg(int index1, int index2) {
	CustomNode* node1 = get(index1);
	CustomNode* node2 = get(index2);
	node1->swapData(node2);
}

//swap(index1,index2)		//Swaps node located at index1 with node at index2
void Slist::swap(int index1, int index2) {
	if(index1 == index2) return;
	else if(index1 > index2) {
		swap(index2, index1);
		return;
	}
	CustomNode* ind1 = get(index1);
	CustomNode* ind2 = get(index2);
	if(ind2 == nullptr) return;
	CustomNode* after2 = get(index2+1);

	CustomNode* before1;
	if(index1 == 0) before1 = head;
	else before1 = get(index1-1);

	if(index1+1 == index2) {
		before1->ptr = ind2;
		ind2->ptr = ind1;
		ind1->ptr = after2;
	} else {
		CustomNode* after1 = get(index1+1);
		CustomNode* before2 = get(index2-1);
		before1->ptr = ind2;
		ind2->ptr = after1;
		before2->ptr = ind1;
		ind1->ptr = after2;
	}
}

// isEmpty()				//Returns true if this list contains no elements.
bool Slist::isEmpty() {
	return head == nullptr;
}

// remove(index)			//Removes the element at the specified index from this list.
void Slist::remove(int index) {
	CustomNode* before;
	if(index == 0) {
		before = head;
	} else {
		before = get(index-1);
	}
	CustomNode* after = get(index+1);
	delete get(index);
	before->ptr = after;
	sizeInt--;
}

// set(index, value)		//Replaces the element at the specified index in this list with a new value.
void Slist::set(int index, CustomNode *value) {
	//FIX
	if(index >= sizeInt) throw 20;
	CustomNode* before = get(index-1);
	value->ptr = before->ptr->ptr;
	delete before->ptr;
	before->ptr = value;
	//insert(index, value);
	//remove(index+1);
}

// size()					//Returns the number of elements in this list.
int Slist::size() {
	return sizeInt;
}

// subList(start, length)	//Returns a new list containing elements from a sub-range of this list.
Slist* Slist::subList(int start, int length) {
	//cout << "checking out of bounds" << endl;
	if(start+length > sizeInt) {
		cout << "Index out of bounds exception: " << start << " " << length  << endl;
		throw 20;
	}
	auto *biglist = new Slist();
	for(int i = start; i<(start+length); ++i) {
		biglist->add(get(i)->copyData());
	}
	//cout << "Sublist: " << biglist->toString() << "\n" << endl;
	return biglist;
}

// toString()				//Converts the list to a printable string representation.
string Slist::toString() {
	string ret;
	int s = 0;
	CustomNode* last = head->ptr;
	while(last->ptr != nullptr) {
		ret += last->toString() + "\n";
		cout << s << " " << last->toString() << endl;
		last = last->ptr;
		s++;
	}
	return ret;
	/*string ret;
	CustomNode* last = head->ptr;
	while(last->ptr != nullptr) {
		ret += last->toString() + "\n";
		last = last->ptr;
	}
	return ret;*/
}

string Slist::toString(CustomNode* ref) {
	int s = 0;
	CustomNode* last = head->ptr;
	while(last->ptr != nullptr) {
		cout << s << " " << last->toString(ref) << endl;
		//cout << last->getCode() << endl;
		last = last->ptr;
		s++;
	}
	return "";
}

void Slist::mapAll(void *fn(CustomNode*)) {
  for(int i = 0; i<sizeInt; ++i)
	fn(get(i));
}

