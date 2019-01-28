#ifndef AIRPORT_SLIST_H
#define AIRPORT_SLIST_H

#include <iostream>
#include <string>
#include <cstring>

using namespace std;

/*

Header file for single linked list class library

*/

class CustomNode {
public:
	CustomNode* ptr;

	virtual string toString() {return "No data (Custom Node type)";}
	virtual void swapData(CustomNode*) {throw;}

	virtual CustomNode* copyData() {throw;}
	virtual string toString(CustomNode*) {return "No data (Custom Node type)";}

	virtual CustomNode *copyDataAndPtr() {
		throw;
	}
};

class Slist {
private:
	int sizeInt;
	CustomNode* tail;
public:

	//Constructor
	Slist();
	//Destructor
	~Slist();

	// add(value)				//Adds a new value to the end of this list.
	void add(CustomNode* value);

	// clear()					//Removes all elements from this list.
	void clear();

	// equals(list)				//Returns true if the two lists contain the same elements in the same order.
	bool equals(Slist* list);

	//get(index)				//Returns the element at the specified index in this list.
	CustomNode* get(int index);

	//insert(index, value)		//Inserts the element into this list before the specified index.
	void insert(int index, CustomNode* value);

	//exchg(index1, index2)		//Switches the payload data of specified indexex.
	void exchg(int index1, int index2);

	//swap(index1,index2)		//Swaps node located at index1 with node at index2
	void swap(int index1, int index2);

	// isEmpty()				//Returns true if this list contains no elements.
	bool isEmpty();

	// mapAll(fn)				//Calls the specified function on each element of the linkedlist in ascending index order.
	void mapAll(void *fn(CustomNode*));

	// remove(index)			//Removes the element at the specified index from this list.
	void remove(int index);

	// set(index, value)		//Replaces the element at the specified index in this list with a new value.
	void set(int index, CustomNode* value);

	// size()					//Returns the number of elements in this list.
	int size();

	// subList(start, length)	//Returns a new list containing elements from a sub-range of this list.
	Slist* subList(int start, int length);

	// toString()				//Converts the list to a printable string representation.
	string toString();

	string toString(CustomNode *ref);

	CustomNode* head;
};


#endif
