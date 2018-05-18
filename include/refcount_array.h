#ifndef _REFCOUNT_ARRAY_H
#define _REFCOUNT_ARRAY_H

// This template is intended for situations where std::vector doesn't work,
// such as when you need access to raw pointers, or when you want to make
// copies without actually duplicating the array.  Memory management is
// handled automatically, but be wary of raw pointers, as they may become
// invalid.

#include <list>		// list<>
#include <sys/types.h>	// size_t
#include <vector>	// vector<>

template<typename _Tp> class _refcount_array {
    public:
	typedef _Tp			value_type;
	typedef value_type *		pointer;
	typedef size_t			size_type;
    private:
	size_type _M_count;	// number of references to this data
	size_type _M_size;	// size of array
	pointer _M_array;
    public:
	explicit _refcount_array(void) : _M_count(1), _M_size(0), _M_array(0) { }
	explicit _refcount_array(const size_type __n) : _M_count(1), _M_size(__n), _M_array() {
		if (_M_size == 0) {
			_M_array = 0;
		} else {
			_M_array = new value_type[_M_size];
		}
	}
	explicit _refcount_array(const std::list<value_type> &__a) : _M_count(1), _M_size(__a.size()), _M_array() {
		if (_M_size == 0) {
			_M_array = 0;
		} else {
			_M_array = new value_type[_M_size];
			typename std::list<value_type>::const_iterator __b(__a.begin());
			const typename std::list<value_type>::const_iterator __end_b(__a.end());
			pointer __c(begin());
			for (; __b != __end_b; ++__b, ++__c) {
				*__c = *__b;
			}
		}
	}
	explicit _refcount_array(const std::vector<value_type> &__a) : _M_count(1), _M_size(__a.size()), _M_array() {
		if (_M_size == 0) {
			_M_array = 0;
		} else {
			_M_array = new value_type[_M_size];
			typename std::vector<value_type>::const_iterator __b(__a.begin());
			const typename std::vector<value_type>::const_iterator __end_b(__a.end());
			pointer __c(begin());
			for (; __b != __end_b; ++__b, ++__c) {
				*__c = *__b;
			}
		}
	}
	~_refcount_array(void) {
		if (_M_array != 0) {
			delete[] _M_array;
		}
	}
	bool empty(void) const {
		return _M_size == 0;
	}
	size_type size(void) const {
		return _M_size;
	}
	void resize(const size_type __n) {
		if (_M_size != __n) {
			_M_size = __n;
			if (_M_array != 0) {
				delete[] _M_array;
			}
			if (_M_size == 0) {
				_M_array = 0;
			} else {
				_M_array = new value_type[_M_size];
			}
		}
	}
	void increment(void) {
		++_M_count;
	}
	bool decrement(void) {
		--_M_count;
		return _M_count == 0;
	}
	pointer array(void) {
		return _M_array;
	}
	pointer begin(void) {
		return _M_array;
	}
	pointer end(void) {
		return _M_array + _M_size;
	}
    private:
	_refcount_array(const _refcount_array<value_type> &);
	_refcount_array<value_type> &operator=(const _refcount_array<value_type> &);
};

template<typename _Tp> class refcount_array {
    public:
	typedef _Tp			value_type;
	typedef value_type *		pointer;
	typedef const value_type *	const_pointer;
	typedef value_type &		reference;
	typedef const value_type &	const_reference;
	typedef size_t			size_type;
    private:
	typedef refcount_array<value_type>		refcount_array_type;
	_refcount_array<value_type> *_M_data;
    public:
	explicit refcount_array(void) : _M_data(new _refcount_array<value_type>) { }
	explicit refcount_array(const size_type __n) : _M_data(new _refcount_array<value_type>(__n)) { }
	explicit refcount_array(const std::list<value_type> &__a) : _M_data(new _refcount_array<value_type>(__a)) { }
	explicit refcount_array(const std::vector<value_type> &__a) : _M_data(new _refcount_array<value_type>(__a)) { }
	refcount_array(const refcount_array_type &__a) : _M_data(__a._M_data) {
		_M_data->increment();
	}
	~refcount_array(void) {
		if (_M_data->decrement()) {
			delete _M_data;
		}
	}
	refcount_array_type &operator=(const refcount_array_type &__a) {
		// make sure we don't already point to this array
		if (_M_data != __a._M_data) {
			// remove reference to old array
			if (_M_data->decrement()) {
				delete _M_data;
			}
			// copy in reference to new array
			_M_data = __a._M_data;
			_M_data->increment();
		}
		return *this;
	}
	refcount_array_type &operator=(const std::list<value_type> &__a) {
		// remove reference to old array
		if (_M_data->decrement()) {
			delete _M_data;
		}
		// copy in new data
		_M_data = new _refcount_array<value_type>(__a);
		return *this;
	}
	refcount_array_type &operator=(const std::vector<value_type> &__a) {
		// remove reference to old array
		if (_M_data->decrement()) {
			delete _M_data;
		}
		// copy in new data
		_M_data = new _refcount_array<value_type>(__a);
		return *this;
	}
	// note that pointers to the array can be invalidated if anything
	// happens to the array, so only keep pointers for short periods
	// where you can guarantee that that won't happen
	pointer array(void) {		 	// direct access to array
		return _M_data->array();
	}
	const_pointer array(void) const {	// direct access to array
		return _M_data->array();
	}
	reference operator[](const size_type __n) {
		return *(_M_data->array() + __n);
	}
	const_reference operator[](const size_type __n) const {
		return *(_M_data->array() + __n);
	}
	reference front(void) {
		return *_M_data->array();
	}
	const_reference front(void) const {
		return *_M_data->array();
	}
	reference back(void) {
		return *(_M_data->array() + _M_data->size() - 1);
	}
	const_reference back(void) const {
		return *(_M_data->array() + _M_data->size() - 1);
	}
	// helper functions that get passed directly to the object
	bool empty(void) const {
		return _M_data->empty();
	}
	void resize(const size_type __n) {
		_M_data->resize(__n);
	}
	size_type size(void) const {
		return _M_data->size();
	}
	pointer begin(void) {
		return _M_data->begin();
	}
	const_pointer begin(void) const {
		return _M_data->begin();
	}
	pointer end(void) {
		return _M_data->end();
	}
	const_pointer end(void) const {
		return _M_data->end();
	}
};

#endif // !_REFCOUNT_ARRAY_H
