#pragma once

template <class Container>
class Node {
public:
	Container* node;
	Node* next;

	Node(Container* _node) {
		node = _node;
		next = nullptr;
	}

	Node() {
		node = nullptr;
		next = nullptr;
	}
};

template <class Container>
class ContainerList {
public:
	Node<Container>* first;
	Node<Container>* last;
	int count;

	ContainerList() : first(nullptr), last(nullptr) { count = 0; }

	bool is_empty() {
		return first == nullptr ? true : false;
	}

	void push_back(const Container &_container) {
		Node<Container> _p = Node<Container>(_container);
		Node<Container>* p = &_p;
		p->node = _container;
		if (is_empty()) {
			first = p;
			last = p;
			count++;
			return;
		}
		last->next = p;
		last = p;
		count++;
	}

	void remove_first() {
		if (is_empty()) return;
		Node<Container>* p = first;
		p->node->remove_data();
		first = p->next;
		delete p;
		count--;
	}

	void remove_last() {
		if (is_empty()) return;
		if (first == last) {
			remove_first();
			count--;
			return;
		}
		Node* p = first;
		while (p->next != last) p = p->next;
		p->next = nullptr;
		last->remove_data();
		delete last;
		last = p;
		count--;
	}

	void clear() {
		if (is_empty()) return;
		while (first != nullptr) {
			remove_first();
		}
		count = 0;
	}

	Container* operator[] (const int index) {
		if (is_empty()) return nullptr;
		Node<Container>* p = first;
		for (int i = 0; i < index; i++) {
			p = p->next;
			if (!p) return nullptr;
		}
		return p->node;
	}
};