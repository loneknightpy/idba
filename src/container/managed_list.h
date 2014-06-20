/**
 * @file managed_list.h
 * @brief ManagedList class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-25
 */

#ifndef __CONTAINER_MANAGED_LIST_H_

#define __CONTAINER_MANAGED_LIST_H_

#include <cstddef>

#include "basic/pool.h"

template <typename T>
struct ManagedListNode
{
    ManagedListNode<T> *next;
    T value;
};

template <typename T>
class ManagedList;

template <typename T>
class ManagedListIterator
{
public:
    typedef T value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef ManagedListNode<T> node_type;
    typedef ManagedList<T> list_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef ManagedListIterator<T> iterator;

    friend class ManagedList<T>;

    ManagedListIterator(node_type *current = NULL)
        : current_(current) {}
    ManagedListIterator(const iterator &iter)
        : current_(iter.current_) {}

    const iterator &operator =(const iterator &iter)
    { current_ = iter.current_; return *this; }

    bool operator ==(const iterator &iter) const
    { return current_ == iter.current_; }
    bool operator !=(const iterator &iter) const
    { return current_ != iter.current_; }

    reference operator *() const { return current_->value; }
    pointer operator ->() const { return &current_->value; }

    const iterator &operator ++()
    { increment(); return *this; }
    iterator operator ++(int)
    { iterator tmp(*this); increment(); return tmp; }

private:
    void increment()
    { if (current_ != NULL) current_ = current_->next; }

    node_type *current_;
};

template <typename T>
class ManagedListConstIterator
{
public:
    typedef T value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef ManagedListNode<T> node_type;
    typedef ManagedList<T> list_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef ManagedListConstIterator<T> const_iterator;

    ManagedListConstIterator(const node_type *current = NULL)
        : current_(current) {}
    ManagedListConstIterator(const const_iterator &iter)
        : current_(iter.current_) {}

    const const_iterator &operator =(const const_iterator &iter)
    { current_ = iter.current_; return *this; }

    bool operator ==(const const_iterator &iter) const
    { return current_ == iter.current_; }
    bool operator !=(const const_iterator &iter) const
    { return current_ != iter.current_; }

    const_reference operator *() const { return current_->value; }
    const_pointer operator ->() const { return &current_->value; }

    const const_iterator &operator ++()
    { increment(); return *this; }
    const_iterator operator ++(int)
    { const_iterator tmp(*this); increment(); return tmp; }

private:
    void increment()
    { if (current_ != NULL) current_ = current_->next; }

    const node_type *current_;
};

/**
 * @brief It is a list impletation which relies on an external Pool instance 
 * for memory manangement.
 *
 * @tparam T
 */
template <typename T>
class ManagedList
{
public:
    typedef T value_type;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;

    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;

    typedef ManagedListNode<T> node_type;
    typedef ManagedList<T> list_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef ManagedListIterator<T> iterator;
    typedef ManagedListConstIterator<T> const_iterator;
    typedef Pool<node_type> node_pool_type;

    ManagedList()
        : size_(0), head_(NULL), pool_(NULL)
    {}
    explicit ManagedList(node_pool_type &pool)
        : size_(0), head_(NULL), pool_(&pool) 
    {}
    ManagedList(const list_type &list)
        : size_(0), head_(NULL), pool_(list.pool_)
    { assign(list); }

    ~ManagedList()
    { clear(); }

    const list_type &operator =(const list_type &list)
    { return assign(list); }

    const list_type &assign(const list_type &list)
    {
        if (this != &list)
        {
            clear();

            node_type *prev = NULL;
            for (node_type *node = list.head_; node; node = node->next)
            {
                node_type *p = pool_->construct();
                p->value = node->value;
                p->next = NULL;

                if (prev == NULL)
                    head_ = p;
                else
                    prev->next = p;
                prev = p;
            }
            size_ = list.size_;
        }

        return *this;
    }

    iterator begin() { return iterator(head_); }
    const_iterator begin() const { return const_iterator(head_); }
    iterator end() { return iterator(); }
    const_iterator end() const { return const_iterator(); }

    void push_front(const value_type &value)
    {
        node_type *p = pool_->construct();
        p->value = value;
        p->next = head_;
        head_ = p;
        ++size_;
    }

    void pop_front()
    {
        node_type *p = head_;
        head_ = head_->next;
        pool_->destroy(p);
        --size_;
    }

    value_type &front() { return head_->value; }
    const value_type &front() const { return head_->value; }

    void remove(const value_type &value)
    {
        node_type *prev = NULL;
        node_type *node = head_;
        while (node)
        {
            if (node->value == value)
            {
                node_type *p = node;
                if (prev == NULL)
                    head_ = node->next;
                else
                    prev->next = node->next;

                node = node->next;
                pool_->destroy(p);
                --size_;
            }
            else
            {
                prev = node;
                node = node->next;
            }
        }
    }

    void erase(iterator p)
    {
        if (p.current_)
        {
            node_type *prev = NULL;
            node_type *node = head_;
            while (node)
            {
                if (node == p.current_)
                {
                    node_type *p = node;
                    if (prev == NULL)
                        head_ = node->next;
                    else
                        prev->next = node->next;

                    node = node->next;
                    pool_->destroy(p);
                    --size_;
                }
                else
                {
                    prev = node;
                    node = node->next;
                }
            }
        }
    }

    iterator find(const value_type &value)
    {
        for (node_type *node = head_; node; node = node->next)
        {
            if (value == node->value)
                return iterator(node);
        }
        return iterator();
    }

    const_iterator find(const value_type &value) const
    {
        for (node_type *node = head_; node; node = node->next)
        {
            if (value == node->value)
                return const_iterator(node);
        }
        return const_iterator();
    }

    node_pool_type &pool() { return *pool_; }
    void set_pool(node_pool_type &pool) { pool_ = &pool; }

    void swap(list_type &list)
    {
        if (this != &list)
        {
            std::swap(size_, list.size_);
            std::swap(head_, list.head_);
            std::swap(pool_, list.pool_);
        }
    }

    size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }

    void clear()
    {
        while (head_)
        {
            node_type *p = head_;
            head_ = head_->next;
            pool_->destroy(p);
        }
        size_ = 0;
    }

private:
    size_type size_;
    node_type *head_;
    node_pool_type *pool_;
};

namespace std
{
template <typename T> inline void swap(ManagedList<T> &x, ManagedList<T> &y)
{ x.swap(y); }
}

#endif

