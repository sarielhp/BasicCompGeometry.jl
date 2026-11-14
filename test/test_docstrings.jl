# Every exported name of the package and its submodules has a docstring. A
# docstring separated from its definition (or written as a bare raw"..."
# string, which Julia does not attach) leaves `?name` empty.
# Base.Docs.hasdoc exists from Julia 1.11; older versions look up the binding in
# the documentation tables of the loaded modules.
function _has_doc(m::Module, n::Symbol)
    isdefined(Base.Docs, :hasdoc) && return Base.Docs.hasdoc(m, n)
    b = Base.Docs.Binding(m, n)
    return any(mod -> haskey(Base.Docs.meta(mod), b), Base.Docs.modules)
end

@testset "Exported names are documented" begin
    M = BasicCompGeometry
    subs = [
        getfield(M, n) for n in names(M; all = true) if isdefined(M, n) &&
            getfield(M, n) isa Module &&
            getfield(M, n) !== M &&
            parentmodule(getfield(M, n)) === M
    ]
    missing_docs = String[]
    for m in (M, subs...), n in names(m)
        n === nameof(m) && continue
        _has_doc(m, n) || push!(missing_docs, "$(nameof(m)).$(n)")
    end
    @test isempty(missing_docs) || (@info "Undocumented exports" missing_docs; false)
end
