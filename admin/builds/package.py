import os

def do_build(context):
    # The environment variables are intended to be set as a build
    # parameter (or otherwise in the Jenkins job configuration).
    version = os.getenv('PACKAGE_VERSION_STRING', 'unknown')
    release = os.getenv('RELEASE', None)
    if not release:
        version += '-dev'

    package_name = 'regressiontests-' + version

    context.make_archive(package_name, use_git=True, prefix=package_name)

    package_name += '.tar.gz'
    package_info = {
            'REGRESSIONTESTS_PACKAGE_FILE_NAME': package_name,
            'REGRESSIONTESTS_PACKAGE_VERSION': version,
            'REGRESSIONTESTS_MD5SUM': context.compute_md5(package_name)
        }
    log_path = context.workspace.get_path_for_logfile('package-info.log')
    context.write_property_file(log_path, package_info)
