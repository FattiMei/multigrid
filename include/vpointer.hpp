#include <stdexcept>


template<typename T>
class PointerVariant {
	public:
		PointerVariant()           : type(Type::Const)  ,   const_ptr(nullptr) {}
		PointerVariant(T *p)       : type(Type::Mutable),   mutable_ptr(p)     {}
		PointerVariant(const T *p) : type(Type::Const)  ,   const_ptr(p)       {}


		const T* get_const_ptr() const {
			return const_ptr;
		}


		T* get_mutable_ptr() const {
			T* result = nullptr;

			switch (type) {
				case Type::Mutable:
					result = mutable_ptr;
					break;

				case Type::Const:
					throw std::runtime_error("Attempt to use a const pointer for non-const purpuoses");
					break;

				// std::unreachable will be available in C++23
			}

			return result;
		}


		PointerVariant<T> operator+(int offset) const {
			PointerVariant<T> result;

			switch (type) {
				case Type::Mutable:
					result = PointerVariant(mutable_ptr + offset);
					break;

				case Type::Const:
					result = PointerVariant(const_ptr + offset);
					break;
			}

			return result;
		}


	private:
		enum class Type {Mutable, Const} type;
		union {
			T* mutable_ptr;
			const T* const_ptr;
		};
};
