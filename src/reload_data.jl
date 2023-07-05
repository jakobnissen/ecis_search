using FASTX
using BioSequences
using BioAlignments

# Load data
function parse_faa(path)
    open(FASTAReader, path) do reader
        map(reader) do record
            (String(description(record)), LongAA(FASTX.sequence(record)))
        end
    end
end

positives_experimental = parse_faa("choices/query_toxins.faa")
positives_mutational = parse_faa("raw/mutational_confirmed.faa")
positives_blast = parse_faa("tmp/possible_toxins.faa")

negatives_experimental = parse_faa("raw/Afp16followingproteins_negative.fasta")
negatives_mutational = parse_faa("raw/mutational_rejected.faa")

negatives_homologous = let
    x = parse_faa("raw/Afp16followingproteins_negative_2.fasta")
    confirmed_identifiers = Set(first.(negatives_experimental))
    filter!(x) do (name, seq)
        !in(name, confirmed_identifiers)
    end
end
negatives_antitoxins = parse_faa("raw/20230313_Antitoxins_2.fasta")

all_positives = vcat(
    positives_experimental,
    positives_mutational,
    positives_blast
)

all_negatives = vcat(
    negatives_experimental,
    negatives_mutational,
    negatives_antitoxins,
    negatives_homologous
)

# Cut off all but leading 20 AA
for pair in (all_positives, all_negatives), (name, seq) in pair
    resize!(seq, 20)
end

# Dereplicate positives
# We don't rerep negatives because we want to be able to discriminate for the small things
MODEL =  AffineGapScoreModel(BioAlignments.BLOSUM62, gap_open = -10, gap_extend = -2)

function dereplicate(v::Vector{<:Tuple{AbstractString, AASeq}}, max_id=0.90)
    remaining = empty(v)
    for (newname, newseq) in v
        push = true
        for (oldname, oldseq) in remaining
            if get_id(oldseq, newseq) > max_id
                push = false
                break
            end
        end
        push && push!(remaining, (newname, newseq))
    end
    remaining
end

function get_id(a, b)
    aas = collect(pairalign(GlobalAlignment(), a, b, MODEL).aln)
    sum(i == j for (i,j) in aas) / length(aas)
end

all_positives_derep = dereplicate(all_positives)

# For negatives, we are less scared of homology, so we just use 100%
# homology reduction
all_negatives_dedup = let
    Dict(v=>k for (k,v) in Dict(v=>k for (k,v) in all_negatives))
end


# Display inter-sequence identity?
for i in 1:length(all_positives_derep)-1, j in i+1:length(all_positives_derep)
    println(get_id(all_positives_derep[i][2], all_positives_derep[j][2]))
end

# Dump the files
for (name, pairs) in [
    ("re_positives.faa", all_positives_derep),
    ("re_negatives.faa", all_negatives_dedup)
]
    open("tmp/"*name, "w") do io
        for (name, seq) in pairs
            println(io, ">", name, '\n', seq)
        end
    end
end 
