# Search blast results
# BLAST'ed conserved proteins against a large bacteria database:
# Identify 50 kbp areas where all 3 proteins are present.

using Intervals
using BioSequences
using FASTX
using BioAlignments

length_of_query = let
    result = Dict{String, Int}()
    for file in ["choices/query_afp.faa", "choices/query_toxins.faa"]
        for record in open(collect, FASTAReader, file)
            result[identifier(record)] = seqsize(record)
        end
    end
    result
end

const BlastHit = @NamedTuple{
    qacc::SubString{String},
    sacc::SubString{String},
    ident::Float64,
    length::Int,
    qstart::Int,
    qend::Int,
    sstart::Int,
    send::Int,
    evalue::Float64,
    bitscore::Float64,
}

function parse_blast(line::String)::BlastHit
    fields = split(line)
    length(fields) == 12 || error("Not exactly 12 fields")
    BlastHit((
        fields[1],
        fields[2],
        parse(Float64, fields[3]) / 100,
        parse(Int, fields[4]; base=10),
        parse(Int, fields[7]; base=10),
        parse(Int, fields[8]; base=10),
        parse(Int, fields[9]; base=10),
        parse(Int, fields[10]; base=10),
        parse(Float64, fields[11]),
        parse(Float64, fields[12]),
    ))
end

read_blast(io::IO)::Vector{BlastHit} = map(parse_blast, eachline(io))

@enum Kind::UInt8 afp1 afp11 afp15 toxin

# These are from the Afp11-specificstrains.fasta and 2 similar FASTA files
KIND_DICT = Dict([
    "AQS36644.1" => afp1,
    "AKH65092.1" => afp1,
    "AUQ43683.1" => afp1,
    "AKH65177.1" => afp1,
    "CAE13960.1" => afp1,
    "AAT48342.1" => afp1,
    "CAQ84339.1" => afp1,
    "AAT48338.1" => afp1,
    "CAE13978.1" => afp1,
    "CAE14001.1" => afp1,
    "KWZ46933.1" => afp1,
    "AKH64624.1" => afp1,
    "CAQ84883.1" => afp1,
    "CAQ85437.1" => afp1,
    "CBW75291.1" => afp1,
    "CAE14019.1" => afp1,
    "AUQ41798.1" => afp1,
    "CAE13997.1" => afp1,
    "CAE14719.1" => afp1,
    "CAE14710.1" => afp1,
    "CAE13982.1" => afp1,
    "AKH64628.1" => afp1,
    "CAE14901.1" => afp1,
    "CAQ84879.1" => afp1,
    "CBW73552.1" => afp1,
    "AWK42277.1" => afp1,
    "CBW76024.1" => afp1,
    "AKH65096.1" => afp1,
    "CBW75294.1" => afp1,
    "CAQ84204.1" => afp1,
    "AKH64611.1" => afp1,
    "CAE13956.1" => afp1,
    "CAE14905.1" => afp1,
    "CAQ84085.1" => afp1,
    "AQS36647.1" => afp1,
    "AKH65181.1" => afp1,
    "CAQ85441.1" => afp1,
    "CAQ84089.1" => afp1,
    "CBW73556.1" => afp1,
    "CAQ84201.1" => afp1,
    "AKH64608.1" => afp1,
    "KWZ46939.1" => afp1,
    "CBW76029.1" => afp1,
    "CAQ84336.1" => afp1,
    "AWK42273.1" => afp1,
    "CAE14023.1" => afp1,
    "CAE13972.1" => afp11,
    "CAQ84330.1" => afp11,
    "KWZ46944.1" => afp11,
    "CAE14732.1" => afp11,
    "CAE13991.1" => afp11,
    "CBW75285.1" => afp11,
    "CAQ84195.1" => afp11,
    "AKH65102.1" => afp11,
    "CAQ84888.1" => afp11,
    "AKH65171.1" => afp11,
    "CAE13950.1" => afp11,
    "AQS36638.1" => afp11,
    "AKH64618.1" => afp11,
    "CAQ85431.1" => afp11,
    "CAQ84095.1" => afp11,
    "CAE14013.1" => afp11,
    "AAT48348.1" => afp11,
    "CBW76018.1" => afp11,
    "AUQ41804.1" => afp11,
    "CAE14895.1" => afp11,
    "AWK42267.1" => afp11,
    "CBW73547.1" => afp11,
    "AKH64602.1" => afp11,
    "CAQ85427.1" => afp15,
    "CAE13986.1" => afp15,
    "AKH65106.1" => afp15,
    "AQS36634.1" => afp15,
    "CAE14009.1" => afp15,
    "AKH65167.1" => afp15,
    "CBW73534.1" => afp15,
    "CAE14741.1" => afp15,
    "CAQ84324.1" => afp15,
    "CAE13968.1" => afp15,
    "CBW75280.1" => afp15,
    "CBW76010.1" => afp15,
    "AWK42264.1" => afp15,
    "CAQ84099.1" => afp15,
    "KWZ46948.1" => afp15,
    "AUQ41807.1" => afp15,
    "AAT48352.1" => afp15,
    "CAQ84192.1" => afp15,
    "CAE14891.1" => afp15,
    "CAQ84892.1" => afp15,
    "AKH64599.1" => afp15,
    "AKH64614.1" => afp15,
    "CAE13946.1" => afp15,
])

minimum_toxin_length::Int = typemax(Int)

toxin_proteins = Dict{String, LongAA}()
for record in open(collect, FASTAReader, "choices/query_toxins.faa")
    KIND_DICT[identifier(record)] = toxin
    global minimum_toxin_length = min(minimum_toxin_length, seqsize(record))
    toxin_proteins[identifier(record)] = FASTX.sequence(LongAA, record)
end

KINDS = unique(values(KIND_DICT))

function filter_hits(hits, min_id, min_cov)
    filter(hits) do hit
        (start, stop) = minmax(hit.qstart, hit.qend)
        cov = (stop - start + 1) / length_of_query[hit.qacc]
        hit.ident ≥ min_id && cov ≥ min_cov
    end
end

afp_blasts = open(read_blast, "tmp/tblastn.tsv")
append!(afp_blasts, open(read_blast, "tmp/tblastn_env_query.tsv"))
afp_blasts = filter_hits(afp_blasts, 0.25, 0.5)

toxin_blasts = open(read_blast, "tmp/tblastn_toxins.tsv")
append!(toxin_blasts, open(read_blast, "tmp/tblastn_env_toxins.tsv"))
toxin_blasts = filter_hits(toxin_blasts, 0.25, 0.5)

blasts = vcat(afp_blasts, toxin_blasts)

# Step one: Reduce search space: Only subjects with all classes
has_all = let
    members = Dict(k => Set{String}() for k in KINDS)
    for row in blasts
        push!(members[KIND_DICT[row.qacc]], row.sacc)
    end
    intersect(values(members)...)
end

# At each subject, create an IntervalSet for each class
class_positions = let
    res = Dict(subject => Dict{Kind, typeof(IntervalSet([1..2]))}() for subject in has_all)
    for row in blasts
        row.sacc in has_all || continue
        (a, b) = minmax(row.send, row.sstart)
        interval = (a-50_000)..(b+50_000)
        d = res[row.sacc]
        kind = KIND_DICT[row.qacc]
        if haskey(d, kind)
            d[kind] = union(d[kind], IntervalSet([interval]))
        else
            d[kind] = IntervalSet([interval])
        end
    end
    res
end

merged_positions = let
    map(collect(class_positions)) do (source, d)
        intervals = collect(values(d))
        result = first(intervals)
        for i in intervals
            result = intersect(result, i)
        end
        (source, result)
    end
end |> Dict
filter!(merged_positions) do (k, v)
    !isempty(v)
end

toxins_inside_intersections = filter(toxin_blasts) do hit
    intervalset = get(merged_positions, hit.sacc, nothing)
    intervalset === nothing && return false
    (a, b) = minmax(hit.send, hit.sstart)
    (a ∈ intervalset) && (b ∈ intervalset)
end

# I downloaded the unique genomes that contain a toxin from NCBI
# to tmp/toxin_hosts.fna

# Run find_genes.py to find genes
if !isfile("tmp/predicted_genes.fna")
    error("You must run `src/find_genes.py` to create tmp/predicted_genes.fna")
end
# Format:
# >MT039196.1_23 # 19747 # 19896 # 1 # ID=1_23;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.493


sequence_of_genome = open(FASTAReader, "tmp/toxin_hosts.fna") do reader
    map(reader) do record
        (identifier(record) => FASTX.sequence(LongDNA{4}, record))
    end |> Dict
end

predicted_genes = let
    # Identifier, [(span, is_reverse_complement)...]
    result = Dict{String, Vector{Tuple{UnitRange, Bool}}}()
    header_regex = r"^([A-Z0-9\.])+_\d+"
    open(FASTAReader, "tmp/predicted_genes.fna") do reader
        for record in reader
            # No genes too small to possibly be a toxin
            seqsize(record) < 0.75 * minimum_toxin_length && continue
            fields = split(description(record), " # ")
            # No partial genes
            occursin(r"partial=00", fields[5]) || continue
            identifier = first(split(first(fields), '_'))
            (start, stop) = (parse(Int, fields[2]), parse(Int, fields[3]))
            @assert start < stop
            element = if fields[4] == "1"
                (start:stop, false)
            else
                # If reverse-complement, we get the opposite strand
                @assert fields[4] == "-1"
                len = length(sequence_of_genome[identifier])
                ((len - start + 1):(len - stop + 1), true)
            end
            push!(get!(valtype(result), result, identifier), element)
        end
    end
    foreach(i -> sort!(unique!(i); by=i -> first(first(i))), values(result))
    result
end

potential_toxins = let
    # Toxin => {Source => (Location, is_reverse_complement)}
    result = Dict{String, Dict{String, Vector{Tuple{UnitRange{Int}, Bool}}}}()
    for hit in toxins_inside_intersections
        toxin = hit.qacc
        source = hit.sacc
        from_to = UnitRange(minmax(hit.sstart, hit.send)...)
        is_rc = hit.sstart == last(from_to)

        # Either the hit begins at the toxin start and so includes the signal
        # peptide, or else we use the predicted genes to find the true toxin 
        # start
        if min(hit.qstart, hit.qend) > 1
            middle_pos = (first(from_to) + last(from_to)) ÷ 2
            i = searchsortedlast(predicted_genes[source], middle_pos; by=i -> first(first(i)))
            if iszero(i)
                continue
            end
            (gene_range, is_rc_hit) = predicted_genes[source][i]
            is_rc == is_rc_hit || continue
            if last(gene_range) < first(from_to)
                # println(predicted_genes[source][i-2:i+2])
                # println(from_to)
                # println(is_rc)
            end
            if length(intersect(gene_range, from_to)) < 0.5 * min(length(gene_range), length(from_to))
                # println(hit)
                # println(gene_range)
                # println(from_to)
                # println()
                continue
            end
        else
            gene_range = from_to
        end
        d = get!(valtype(result), result, toxin)
        push!(get!(valtype(d), d, source), (gene_range, is_rc))
    end
    result
end

MODEL = AffineGapScoreModel(BioAlignments.BLOSUM62, gap_open = -25, gap_extend = -2)

# Sanity check
check = []
dbg = Ref{Any}(nothing)
for (toxin, d) in potential_toxins
    for (source, ranges) in d
        for (range, is_rc) in ranges
            dna = sequence_of_genome[source][range]
            if is_rc
                dna = reverse_complement(dna)
            end
            aa_toxin = toxin_proteins[toxin]
            aa = translate(dna)
            if last(aa) == AA_Term
                aa = aa[1:end-1]
                range = first(range):last(range)-3
            end
            if length(findall(AA_Term, aa)) > 0
                error()
            end
            aas = collect(pairalign(GlobalAlignment(), aa, aa_toxin, MODEL).aln)
            id = sum(i == j for (i,j) in aas) / length(aas)
            id < 0.25 && continue
            push!(check, (toxin, source, range, id, is_rc))
        end
    end
end

# Dump to file
open(FASTAWriter, "tmp/possible_toxins.fna") do fna
    open(FASTAWriter, "tmp/possible_toxins.faa") do faa
        for (toxin, source, span, id, is_rc) in check
            dna = sequence_of_genome[source][span]
            if is_rc
                dna = reverse_complement(dna)
            end
            aa = translate(dna)
            if last(aa) == AA_Term
                aa = aa[1:end-1]
            end
            if length(findall(AA_Term, aa)) > 0
                error()
            end
            id = "$(toxin)_$(source)_$(span)_id=$(round(100*id; digits=1))"
            write(fna, FASTARecord(id, dna))
            write(faa, FASTARecord(id, aa))
        end
    end
end
