
#include <vector>
#include <array>
#include <iterator>

#include <Eigen/Geometry>

template<class T, int DIM, class A = std::allocator<T>>
class Tensor
{
public:
	static const int DIMENSIONS = DIM;
	typedef size_t Index;
//	typedef std::array<size_t, DIM> Position;
//	typedef std::array<size_t, DIM> Size;
	typedef Eigen::Array<int,DIM,1> Position;
	typedef Eigen::Array<int,DIM,1> Size;
	typedef A Allocator;
	typedef std::vector<T, A> Data;
	typedef T Value;
	typedef T& Reference;
	typedef const T& ConstReference;
	typedef typename Data::iterator Iterator;
	typedef typename Data::const_iterator ConstIterator;

	class Subblock
	{
	public:
		typedef std::pair <Position, Size> Slice;

		template<class B>
		class IteratorBase : public std::iterator < std::bidirectional_iterator_tag, T, std::ptrdiff_t, T*, T& >
		{
		public:
			IteratorBase() :
				_block(nullptr),
				_index(0)
			{}

			IteratorBase(B block, Index index) :
				_block(block),
				_index(index)
			{}

			IteratorBase(const IteratorBase&) = default;
			IteratorBase &operator=(const IteratorBase&) = default;

			bool operator==(const IteratorBase &other) const { return _block == other._block && _index == other._index; }
			bool operator!=(const IteratorBase &other) const { return !(*this == other); }

			T* operator->() { return _block ? &((*_block)[_index]) : nullptr; }
			const T* operator->() const { return _block ? &((*_block)[_index]) : nullptr; }

			T& operator*() { return (*_block)[_index]; }
			const T& operator*() const { return (*_block)[_index]; }

			IteratorBase &operator ++ ()
			{
			    ++_index;
				return *this;
			}
			IteratorBase operator ++ (int) { IteratorBase tmp(*this); ++(*this); return tmp; }
			IteratorBase &operator -- ()
			{
			    --_index;
				return *this;
			}
			IteratorBase operator -- (int) { IteratorBase tmp(*this); ++(*this); return tmp; }

			void swap(IteratorBase &other) { IteratorBase tmp(*this); *this = other; other = tmp; }

		private:
			B _block;
			Index _index;
		};

		typedef IteratorBase<const Subblock*> ConstIterator;
		typedef IteratorBase<Subblock*> Iterator;

		Subblock(Tensor &tensor, const Position &position, const Size &size) :
			_tensor(&tensor),
			_position(position),
			_size(size)
		{}

		Subblock(Tensor &tensor, const Slice &slice) :
			Subblock(tensor, std::get<0>(slice), std::get<1>(slice))
		{}

		Subblock(const Subblock&) = default;
		Subblock(Subblock&&) = default;

		Subblock &operator=(const Subblock&) = default;
		Subblock &operator=(Subblock&&) = default;

		Position positionAtIndex(Index index) const { return Tensor::indexToPosition(index, _size); }
		Index indexAtPosition(const Position &position) const { return Tensor::positionToIndex(position, _size); }

		Reference      operator[](Index index)       { assert(isValid()); return (*_tensor)[absoluteIndex(index)]; }
		ConstReference operator[](Index index) const { assert(isValid()); return (*_tensor)[absoluteIndex(index)]; }

		Reference      operator[](const Position &position)       { assert(isValid()); return (*_tensor)[absolutePosition(position)]; }
		ConstReference operator[](const Position &position) const { assert(isValid()); return (*_tensor)[absolutePosition(position)]; }

		const Position &position() const { return _position; }
		const Size &size() const { return _size; }
		Index numValues() const { return Tensor::numValues(_size); }

		bool isValid() const
		{
			if (numValues() == 0) return false;

            if ((_position < 0).any() || (_position >= _size).any())
                return false;

//			const auto &tensorSize = _tensor->size();
//			for (int i = 0; i < _position.size(); ++i)
//				if (_position[i] + _size[i] - 1 >= tensorSize[i])
//					return false;

			return true;
		}

		Iterator begin() { return Iterator(this, 0); }
		ConstIterator begin() const { return Iterator(this, 0); }
		ConstIterator cbegin() const { return begin(); }

		Iterator end() { return Iterator(this, numValues()); }
		ConstIterator end() const { return Iterator(this, numValues()); }
		ConstIterator cend() const { return end(); }

	private:
		void toAbsolutePosition(Position &position) const
		{
			for (int i = 0; i < position.size(); ++i)
				position[i] += _position[i];
		}

		Position absolutePosition(const Position &position) const
		{
			Position result = position;
			toAbsolutePosition(result);
			return result;
		}

		Index absoluteIndex(Index index) const
		{
			Position position = positionAtIndex(index);
			toAbsolutePosition(position);
			index = Tensor::positionToIndex(position, _tensor->size());
			return index;
		}

		Tensor *_tensor;
		Position _position;
		Size _size;
	};

	Tensor(Allocator allocator = Allocator()):
		_size(),
		_data(allocator)
	{
	    clear();
	}

	Tensor(const Size &size, const T &initialValue = T(), Allocator allocator = Allocator()) :
		_size(size),
		_data(allocator)
	{
		_data.assign(numValues(), initialValue);
	}

	Tensor(const Tensor&) = default;
	Tensor(Tensor &&) = default;

	Tensor &operator=(const Tensor&) = default;
	Tensor &operator=(Tensor &&) = default;

    inline bool isEmpty() const { return numValues() == 0; }
	inline static Index numValues(const Size &size) noexcept { return size.prod(); }
	inline Index numValues() const noexcept { return numValues(_size); }
	inline const Size &size() const noexcept { return _size; }

	void resize(const Size &size, const T &initialValue = T())
	{
		Data newData(numValues(size), initialValue, _data.get_allocator());

		Size intersection(size);
		for (int i = 0; i < intersection.size(); ++i)
			intersection[i] = std::min(intersection[i], _size[i]);

		int numCellsInIntersection = numValues(intersection);
		for (int i = 0; i < numCellsInIntersection; ++i)
		{
			auto position = indexToPosition(i, intersection);
			auto oldIndex = positionToIndex(position, _size);
			auto newIndex = positionToIndex(position, size);
			newData[newIndex] = std::move(_data[oldIndex]);
		}
		_size = size;
		_data = std::move(newData);
	}

	void clear() noexcept
	{
	    for (int i = 0; i < _size.size(); ++i)
            _size[i] = 0;

		_data.clear();
	}

	void assign(const Size &size, const T &value = T())
	{
		_data.assign(numValues(size), value);
		_size = size;
	}

	template<class Iterator>
	void assign(const Size &size, Iterator begin, Iterator end)
	{
		int bufferSize = numValues(size);
		if (std::distance(begin, end) < bufferSize)
			throw std::invalid_argument("Provided iterators don't yield enough results to fill the tensor.");

		_size = size;
		_data.assign(begin, begin + numValues());
	}

	Reference      at(size_t index)       { assertCorrectIndex(index); return _data.at(index); }
	ConstReference at(size_t index) const { assertCorrectIndex(index); return _data.at(index); }

	Reference      at(const Position &position)       { assertCorrectPosition(position); return _data.at(indexAtPosition(position)); }
	ConstReference at(const Position &position) const { assertCorrectPosition(position); return _data.at(indexAtPosition(position)); }

	Reference      operator[](size_t index)       { return _data[index]; }
	ConstReference operator[](size_t index) const { return _data[index]; }

	Reference      operator[](const Position &position)       { return _data[indexAtPosition(position)]; }
	ConstReference operator[](const Position &position) const { return _data[indexAtPosition(position)]; }

	Iterator      begin()        { return _data.begin(); }
	Iterator      end()          { return _data.end(); }
	ConstIterator begin()  const { return _data.begin(); }
	ConstIterator end()    const { return _data.end(); }
	ConstIterator cbegin() const { return _data.cbegin(); }
	ConstIterator cend()   const { return _data.cend(); }

	Iterator      rbegin()        { return _data.rbegin(); }
	Iterator      rend()          { return _data.rend(); }
	ConstIterator rbegin()  const { return _data.rbegin(); }
	ConstIterator rend()    const { return _data.rend(); }
	ConstIterator rcbegin() const { return _data.crbegin(); }
	ConstIterator rcend()   const { return _data.crend(); }

	Data &data() noexcept { return _data; }
	const Data &data() const noexcept { return _data; }

	static Position indexToPosition(Index index, const Size &size)
	{
        assert(size.prod() > 0);

		Position position;

		size_t modulo = 1;
		size_t divisor = 1;

		//Example: 3D-Case:
		//pos.x = index % size.x;
		//pos.y = index % (size.x * size.y) / size.x;
		//pos.z = index / (size.x * size.y);

		for (int i = 0; i < size.size(); ++i)
		{
			const size_t &s = size[i];
			modulo *= s;
			position[i] = (index % modulo) / divisor;
			divisor *= size[i];
		}
		return position;
	}
	Position positionAtIndex(Index index) const { return indexToPosition(index,_size); }

	static Index positionToIndex(const Position &position, const Size &size)
	{
	    assertCorrectPosition(position,size);

		Index index = 0;
		size_t factor = 1;

		//Example: 3D-Case:
		//pos.x + pos.y * size.x + pos.z * size.x * size.y;

		for (int i = 0; i < size.size(); ++i)
		{
			index += position[i] * factor;
			factor *= size[i];
		}
		return index;
	}
	Index indexAtPosition(const Position &position) const { return positionToIndex(position, _size); }

	bool operator==(const Tensor &other) const { return _size == other._size && _data == other._data; }
	bool operator!=(const Tensor &other) const { return !(*this == other); }

private:
    inline static void assertCorrectPosition(const Position &position, const Size &size)
    {
        assert((position >= 0).all());
        assert((position < size).all());
    }
    inline void assertCorrectPosition(const Position &position) { assertCorrectPosition(position,_size); }

    inline static void assertCorrectIndex(Index index, const Size &size)
    {
        assert(index >= 0);
        assert(index <= size.prod());
    }
    inline void assertCorrectIndex(Index index) { assertCorrectIndex(index,_size); }

	Size _size;
	Data _data;
};
