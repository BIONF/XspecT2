:::mermaid
classDiagram
    ProbabilisticFilterModel <|-- ProbabilisticFilterSVMModel
    ProbabilisticFilterModel <|-- ProbabilisticSingleFilterModel
    ProbabilisticFilterModel : String filter_display_name
    ProbabilisticFilterModel : String author
    ProbabilisticFilterModel : String author_email
    ProbabilisticFilterModel : String filter_type
    ProbabilisticFilterModel : Path base_path
    ProbabilisticFilterModel : dict display_names
    ProbabilisticFilterModel : int k
    ProbabilisticFilterModel : float fpr
    ProbabilisticFilterModel : int num_hashes
    ProbabilisticFilterModel : dict filters
    ProbabilisticFilterModel : dict non_distinct_kmer_counts
    ProbabilisticFilterModel : __init__(k, filter_display_name, author, author_email, filter_type, base_path, fpr, num_hashes)
    ProbabilisticFilterModel : __dict__()
    ProbabilisticFilterModel : to_dict()
    ProbabilisticFilterModel : slug()
    ProbabilisticFilterModel : fit(dir_path, display_names)
    ProbabilisticFilterModel: calculate_hits(sequence, filter_ids)
    ProbabilisticFilterModel : predict(sequence_input, filter_ids)
    ProbabilisticFilterModel : filter(sequences, threshold, filter_ids)
    ProbabilisticFilterModel : save()
    ProbabilisticFilterModel : load(path)
    ProbabilisticFilterModel : _convert_cobs_result_to_dict(cobs_result)
    ProbabilisticFilterModel : _count_kmers(sequence_input)
    ProbabilisticFilterModel : _get_cobs_index_path()

    class ProbabilisticFilterSVMModel {
        SVM svm
        String svm_kernel
        float svm_c
        __init__(..., kernel, c)
        to_dict()
        set_svm_params(kernel, c)
        fit(dir_path, svm_path, display_names)
        predict(sequence_input, filter_ids) return_by_display_names)
        load(path)
        _get_svm(id_keys)

    }

    class ProbabilisticSingleFilterModel {
        Bloom filter
        __init__(...)
        fit(file_path, display_name)
        calculate_hits(sequence)
        load(path)
        _generate_kmers(sequence)
    }
:::