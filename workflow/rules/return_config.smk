
rule return_config:
    output:
        ws_path("config/config_used.yaml")
    params:
        configfile_path=config_copy
    run:
        import shutil
        shutil.copyfile(params.configfile_path, output[0])
