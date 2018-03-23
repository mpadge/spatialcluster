/* binary search tree adapted from
 * http://www.bogotobogo.com/cplusplus/binarytree.php,
 * and used here just to extract minimal and maximal values 
 * Note that although this might be more neatly done by embedding the functions
 * within the struct/class def, it works by recursively adding pointers to new
 * instances of same struct/class. Because it will be ultimately embedded within
 * a meta-struct object, it is therefore actually easier not to embed the
 * function defs, at the price of a little untidiness.
 * */

#include <iostream>

template <typename T>
struct Tree
{
	T data;
	Tree *left;
	Tree *right;  
	Tree *parent;  
};

template <typename T>
Tree <T> *treeNewNode (T data) 
{
	Tree <T> *node = new Tree <T>;
	node->data = data;
	node->left = nullptr;
	node->right = nullptr;
	node->parent = nullptr;

	return node;
}

template <typename T>
void treeInsertNode (Tree <T> *node, T dat)
{
	Tree <T> *newNode = new Tree <T>;
	newNode->data = dat;
	newNode->left = nullptr;
	newNode->right = nullptr;
	while (node)
    {
	    if (dat <= node->data)
        {
		    if (node->left == nullptr)
            {
			    newNode->parent = node;
		        node->left = newNode;
		        return;
		    }
		    else
                node = node->left;
        }
	    else
        {
		    if (node->right == nullptr)
            {
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
Tree <T> *treeGetNode (Tree <T> *node, T dat)
{
	if (node == nullptr)
        return node;

	if (dat < node->data) 
		return treeGetNode (node->left,dat);
	else if ( dat > node->data)
		return treeGetNode (node->right, dat);
	else
		return node;
}

template <typename T>
int treeSize (struct Tree <T> *node)
{
    if (node == nullptr)
        return 0;
    else
        return treeSize (node->left) + 1 + treeSize (node->right);
}

template <typename T>
T treeMin (Tree <T> *node)
{
	while (node->left)
		node = node->left;

	return node->data;
}

template <typename T>
T treeMax (Tree <T> *node)
{
	while (node->right)
		node = node->right;

	return node->data;
}

template <typename T>
Tree <T> *treeSuccesorInOrder (Tree <T> *node)
{
    /* if the node has right child, seccessor is Tree-Minimum */
    if (node->right != nullptr)
        return node->right;

    Tree <T> *y = node->parent;
    while (y != nullptr && node == y->right)
    {
        node = y;
        y = y->parent;
    }
    return y;
}


template <typename T>
T treePredecessorInOrder (Tree <T> *node)
{
	if (node->left) 
		return treeMax (node->left);

	Tree <T> *y = node->parent;
	while (node == y->left)
    {
		node = y;
		y = y->parent;
	}

	return y->data;
}

template <typename T>
void treeDeleteNode (Tree <T> *root, T dat)
{
	Tree <T> *node = nullptr, *p = nullptr,
         *child = nullptr, *pred = nullptr;
	node = treeGetNode (root, dat);

	if (node->left == nullptr && node->right == nullptr)
    {
		if (node->parent) p = node->parent;
		if (node == p->left) 
			p->left = nullptr;
		else
			p->right = nullptr;
		delete node;
		return;
	}

	if (node->left && node->right)
    {
		T ch_pred = treePredecessorInOrder (node);
		pred = treeGetNode (root, ch_pred);
		if (pred->parent->left == pred) 
			pred->parent->left = nullptr;
		else if (pred->parent->right == pred) 
			pred->parent->right = nullptr;
		node->data = pred->data;
		delete pred;
		return;
	}

	if (node->left) 
		child = node->left;
	else if (node->right)
		child = node->right;
	p = node->parent;
	if (p->left && p->left == node)
    {
		p->left = child;
	}
	else if (p->right && p->right == node)
    {
		p->right = child;
	}
	child->parent = p;
	delete node;
}

template <typename T>
void treeClear (Tree <T> *node)
{
    if (node != nullptr)
    {
        treeClear (node->left);
        treeClear (node->right);
        delete node;
    }
}

/* affirm no memory leaks:
#include <random>

// clang++ -fsanitize=undefined bst.cpp -o junk
// valgrind --tool=memcheck --leak-check=full ./junk

int main()
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

	//Tree <double> *root = treeNewNode(1.0);
	Tree <double> *tree;
    
    for (int i = 0; i < 1000; i++)
    {
        if (i == 0)
            tree = treeNewNode (distribution (generator));
        else
            insertTreeNode(tree, distribution (generator));
    }

    std::cout << "tree (min, max) = (" << treeMin (tree) << ", " <<
        treeMax (tree) << ")" << std::endl;
    
    treeDeleteNode (tree, treeMin (tree));
    treeDeleteNode (tree, treeMax (tree));
    std::cout << "  ---> (" << treeMin (tree) << ", " <<
        treeMax (tree) << ")" << std::endl;

    treeClear(tree);

	return 0;
}
*/
