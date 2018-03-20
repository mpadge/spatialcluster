#pragma once

#include <map>
#include <algorithm> // min_element
#include <iostream>
#include <deque>
#include <climits>
#include <vector>
#include <random>
#include <ctime>


/* Binary Tree */

template <typename T>
struct Tree
{
	T data;
	Tree *left = nullptr;
	Tree *right = nullptr;  
	Tree *parent = nullptr;  
};

template <typename T>
Tree <T> *newTreeNode(T data) 
{
	Tree <T> *node = new Tree <T>;
	node->data = data;
	node->left = NULL;
	node->right = NULL;
	node->parent = NULL;

	return node;
}

template <typename T>
void insertTreeNode(Tree <T> *node, T data)
{
	Tree <T> *newNode = new Tree <T>;
	newNode->data = data;
	newNode->left = NULL;
	newNode->right = NULL;
	while(node) {
	    if(data <= node->data) {
		    if(node->left == NULL) {
			    newNode->parent = node;
		        node->left = newNode;
		        return;
		    }
		    else
		    node = node->left;
	    }
	    else {
		    if(node->right == NULL) {
				newNode->parent = node;
		        node->right = newNode;
		        return;
		    }
		    else
		        node = node->right;
	    }
	}
}

template <typename T>
Tree <T> *getNode(Tree <T> *node, T data)
{
	if(node == NULL) return node;
	if(data < node->data) 
		return getNode(node->left,data);
	else if( data > node->data)
		return getNode(node->right, data);
	else
		return node;
}

template <typename T>
T treeMin(Tree <T> *node)
{
	while(node->left) {
		node = node->left;
	}
	return node->data;
}

template <typename T>
T treeMax(Tree <T> *node)
{
	while(node->right) {
		node = node->right;
	}
	return node->data;
}

template <typename T>
void printTreeInOrder(Tree <T> *node)
{
	if(node == NULL) return;

	printTreeInOrder(node->left);
    std::cout << node->data << " ";
	printTreeInOrder(node->right);
}

template <typename T>
void deleteKey(Tree <T> *root, T data)
{
	Tree <T> *node, *p, *child, *pred;
	// get the node for the key
	node = getNode(root, data);

	// leaf node - just delete the node
	if(node->left == NULL && node->right == NULL) {
		if(node->parent) p = node->parent;
		if(node == p->left) 
			p->left = NULL;
		else
			p->right = NULL;
		delete node;
		return;
	}

	// two children - replace it with its predecessor and delete
	if(node->left && node->right) {
		T ch_pred = predecessorInOrder(node);
		pred = getNode(root, ch_pred);
		if(pred->parent->left == pred) 
			pred->parent->left = NULL;
		else if(pred->parent->right == pred) 
			pred->parent->right = NULL;
		node->data = pred->data;
		delete pred;
		return;
	}

	// only one child 
	// replace it with its child and delete
	if(node->left) 
		child = node->left;
	else if(node->right)
		child = node->right;
	p = node->parent;
	if(p->left && p->left == node) {
		p->left = child;
	}
	else if (p->right && p->right == node) {
		p->right = child;
	}
	child->parent = p;
	delete node;
	return;
}

template <typename T>
void clear(Tree <T> *node)
{
    if(node != NULL) {
        clear(node->left);
        clear(node->right);
        delete node;
    }
}

