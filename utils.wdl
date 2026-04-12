version 1.0

task get_julia_cmd {
    input {
        String use_sysimage = "true"
        String threads = "auto"
    }
    command <<<
        julia_cmd_string="julia --project=/opt/decodeme --startup-file=no"
        if [[ "~{use_sysimage}" == "true" ]]; then
            julia_cmd_string+=" --sysimage=/opt/decodeme/sysimage.so"
        fi
        if [[ "~{threads}" == "auto" ]]; then
            julia_cmd_string+=" --threads=auto"
        fi
        echo "$julia_cmd_string"
    >>>

    output {
        String julia_cmd = read_string(stdout())
    }
}