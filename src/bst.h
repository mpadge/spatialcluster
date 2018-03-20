/* binary search tree adapted from
 * http://www.bogotobogo.com/cplusplus/binarytree.php,
 * and used here just to extract minimal and maximal values */

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
Tree <T> *newTreeNode (T data) 
{
	Tree <T> *node = new Tree <T>;
	node->data = data;
	node->left = nullptr;
	node->right = nullptr;
	node->parent = nullptr;

	return node;
}

template <typename T>
void insertTreeNode (Tree <T> *node, T dat)
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
Tree <T> *getNode (Tree <T> *node, T dat)
{
	if (node == nullptr)
        return node;

	if (dat < node->data) 
		return getNode (node->left,dat);
	else if ( dat > node->data)
		return getNode (node->right, dat);
	else
		return node;
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
T predecessorInOrder (Tree <T> *node)
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
void deleteKey (Tree <T> *root, T dat)
{
	Tree <T> *node, *p, *child, *pred;
	node = getNode (root, dat);

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
		T ch_pred = predecessorInOrder (node);
		pred = getNode (root, ch_pred);
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
void clear (Tree <T> *node)
{
    if (node != nullptr)
    {
        clear (node->left);
        clear (node->right);
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

	//Tree <double> *root = newTreeNode(1.0);
	Tree <double> *tree;
    
    for (int i = 0; i < 1000; i++)
    {
        if (i == 0)
            tree = newTreeNode (distribution (generator));
        else
            insertTreeNode(tree, distribution (generator));
    }

    std::cout << "tree (min, max) = (" << treeMin (tree) << ", " <<
        treeMax (tree) << ")" << std::endl;
    
    deleteKey (tree, treeMin (tree));
    deleteKey (tree, treeMax (tree));
    std::cout << "  ---> (" << treeMin (tree) << ", " <<
        treeMax (tree) << ")" << std::endl;

    clear(tree);

	return 0;
}
*/
