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
    log_path = context.workspace.get_path_for_logfile('package-info.log')
    context.write_package_info(log_path, Project.REGRESSIONTESTS, package_name, version)
