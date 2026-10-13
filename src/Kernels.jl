"""
    AbstractKernel

Arithmetic policy for geometric predicates. Ordinary calls use
`FloatingKernel()`. Algorithms accepting `kernel=` must forward that same
object to every predicate; they must not reconstruct a default internally.

Custom kernels implement the kernel-first predicate methods they support,
for example `turn_sign(kernel, p, q, r)` and `compare_xy(kernel, p, q)`.
Unsupported methods fail rather than silently using floating-point arithmetic.
Kernels select predicate arithmetic, not the storage type of geometric objects.
Persistent structures should store the concrete kernel type with their state;
cached predicates from different kernels must not be mixed.
"""
abstract type AbstractKernel end

"""
    FloatingKernel()
    FloatingKernel(Float32)

Ordinary floating-point predicate arithmetic, using `Float64` by default.
Input coordinates are converted before arithmetic; integer overflow in the
original coordinate type is therefore avoided, but rounding and floating-point
overflow can still affect predicate results. This is not a certified filter.
"""
struct FloatingKernel{T<:AbstractFloat} <: AbstractKernel end
FloatingKernel() = FloatingKernel{Float64}()
FloatingKernel(::Type{T}) where {T<:AbstractFloat} = FloatingKernel{T}()

"""
    ExactKernel()

Exact predicates on integer, rational, and finite binary floating-point inputs.
Integers are widened to `BigInt` before arithmetic. Rational and floating-point
inputs use `Rational{BigInt}`. A floating-point input means its stored binary
value, not the decimal value that may have been intended by the caller.

This kernel does not provide exact geometric constructions: an intersection
previously rounded to floating point remains a rounded input. No epsilon or
fixed-precision `BigFloat` calculation is used to decide signs.
"""
struct ExactKernel <: AbstractKernel end

const BuiltinKernel = Union{FloatingKernel,ExactKernel}

@inline function _kernel_scalar(::FloatingKernel{T}, x::Real) where {T}
    value = T(x)
    isfinite(value) || throw(DomainError(x, "predicate inputs must be finite"))
    return value
end

_kernel_scalar(::ExactKernel, x::Integer) = BigInt(x)
function _kernel_scalar(::ExactKernel, x::Rational)
    isfinite(x) || throw(DomainError(x, "predicate inputs must be finite"))
    return Rational{BigInt}(x)
end
function _kernel_scalar(::ExactKernel, x::AbstractFloat)
    isfinite(x) || throw(DomainError(x, "predicate inputs must be finite"))
    return Rational{BigInt}(x)
end

@inline _predicate_sign(value) = Int(sign(value))

export AbstractKernel, FloatingKernel, ExactKernel
