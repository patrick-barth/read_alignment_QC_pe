process final{
    output:
    path("metadata.txt")

    script:
    """
    cat <<EOF > metadata.txt
    Author: ${params.manifest.author}
    Pipeline version: ${params.manifest.version}
    Working directory: ${workflow.workDir}
    User name: ${workflow.userName}
    Command line: ${workflow.commandLine}
    Container engine: ${workflow.containerEngine}    
    Containers used: ${workflow.container}
    Git repository: ${workflow.repository}
    Repository revision: ${workflow.revision}
    """
}
