#include <stdexcept>


template<typename T>
class PointerVariant {
	public:
		PointerVariant(T *p)       : type(Type::Mutable), mutable_ptr(p) {}
		PointerVariant(const T *p) : type(Type::Const)  ,   const_ptr(p) {}


		const T* get_const_ptr() const {
			return const_ptr;
		}


		T* get_mutable_ptr() const {
			switch (type) {
				case Type::Mutable:
					return mutable_ptr;

				case Type::Const:
					throw std::runtime_error("Attempt to use a const pointer for non-const purpuoses");
					return nullptr;

				// std::unreachable will be available in C++23
			}
		}


		PointerVariant operator+(int offset) const {
			switch (type) {
				case Type::Mutable:
					return PointerVariant(mutable_ptr + offset);

				case Type::Const:
					return PointerVariant(const_ptr + offset);
			}
		}


	private:
		const enum class Type {Mutable, Const} type;
		union {
			T* mutable_ptr;
			const T* const_ptr;
		};
};
