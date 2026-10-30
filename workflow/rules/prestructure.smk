#=============================================================================
# MODULE 2: POPULATION STRUCTURE
#=============================================================================

rule snmf:
    """Run sNMF across K range."""
    input:  geno = W['geno']
    output: W['snmf']
    params: ks = K_START, ke = K_END, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}prestructure/snmf.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/snmf.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule pca_plot:
    """Generate PCA and Tracy-Widom plots (also creates LEA PCA projections/eigenvalues)."""
    input:  lfmm = W['lfmm'], meta = O['metadata']
    output:
        pca = O['pca'],
        pca_svg = O['pca_svg'],
        tracy = O['tracy'],
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    params: inter_dir = INTER
    log:    f"{LOGDIR}prestructure/pca_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca.R \
            {input.lfmm} {input.meta} {output.pca} {output.tracy} {params.inter_dir} > {log} 2>&1
        """

rule cross_entropy_plot:
    """Plot cross-entropy for K selection."""
    input:  snmf = W['snmf']
    output: O['cross_entropy']
    params: ks = K_START, ke = K_END, inter_dir = INTER
    log:    f"{LOGDIR}prestructure/cross_entropy.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_cross_entropy.R \
            {input.snmf} {params.ks} {params.ke} \
            {output} {params.inter_dir} > {log} 2>&1
        """

rule tracy_widom_table:
    """Tracy-Widom test of every PCA axis as a table (6bf499), with LEA's own routine on the
    eigenvalues pca_plot wrote. Its own rule: pca_plot's projections feed EMMAX, so a new output
    there would rerun every scan on an existing tree."""
    input:  eigenvalues = W['pca_eigenvalues']
    output: O['tracy_widom_table']
    log:    f"{LOGDIR}prestructure/tracy_widom_table.log"
    shell:
        """
        Rscript /pipeline/scripts/tracy_widom_table.R \
            {input.eigenvalues} {output} > {log} 2>&1
        """

rule cross_entropy_table:
    """sNMF cross-entropy of every K and run as a table (f59f46). Its own rule: the
    cross-entropy plot is an input of write_summary, so a new output on cross_entropy_plot would
    rerun both on an existing tree."""
    input:  snmf = W['snmf']
    output: O['cross_entropy_table']
    params: ks = K_START, ke = K_END
    log:    f"{LOGDIR}prestructure/cross_entropy_table.log"
    shell:
        """
        Rscript /pipeline/scripts/cross_entropy_table.R \
            {input.snmf} {params.ks} {params.ke} {output} > {log} 2>&1
        """

rule extract_clusters:
    """Extract Q-matrix (cluster assignments) for specific K."""
    input:  snmf = W['snmf'], meta = O['metadata']
    output: clusters_table("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k
    log:    f"{LOGDIR}prestructure/extract_clusters_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/extract_clusters.R \
            {input.snmf} {input.meta} {params.k} {output} > {log} 2>&1
        """

rule structure_barplot:
    """Generate structure barplot for specific K."""
    input:  clusters = clusters_table("{k}")
    output: structure_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, inter_dir = INTER
    log:    f"{LOGDIR}prestructure/structure_plot_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_structure.R \
            {input.clusters} {params.k} {output} {params.inter_dir} > {log} 2>&1
        """

rule pca_structure_plot:
    """Generate PCA with structure pie charts for specific K (uses LEA PCA)."""
    input:
        clusters = clusters_table("{k}"),
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    output: pca_struct_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, inter_dir = INTER
    log:    f"{LOGDIR}prestructure/pca_structure_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca_structure.R \
            {input.clusters} {input.projections} {input.eigenvalues} \
            {params.k} {output} {params.inter_dir} > {log} 2>&1
        """

rule pop_diff_test:
    """Population differentiation test for specific K."""
    input:  snmf = W['snmf']
    output: pop_diff_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params:
        k = lambda wc: wc.k,
        ploidy = PLOIDY,
        inter_dir = INTER
    log:    f"{LOGDIR}prestructure/pop_diff_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pop_diff.R \
            {input.snmf} {params.k} {params.ploidy} \
            {output} {params.inter_dir} > {log} 2>&1
        """
