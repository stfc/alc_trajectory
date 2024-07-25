pipeline {
    agent any

    stages {
        stage('Checkout'){
            steps {
                checkout scm
            }
        }    
        stage('Build and Test GCC-Debug') {
            steps {
                sh '''
                module load cmake
                folder="test-gcc-debug"
                rm -rf $folder && mkdir $folder && pushd $folder
                FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON
                make
                ctest --output-on-failure
                '''
            }
        } 
        stage('Build and Test GCC-Release') {
            steps {
                sh '''
                module load cmake
                folder='test-gcc'
                rm -rf $folder && mkdir $folder && pushd $folder
                FC=gfortran cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON
                make
                ctest --output-on-failure
                '''
            }
        }
        stage('Build and Test Intel-Debug') {
            steps {
                sh '''
                module load cmake
                module load intel
                folder='test-ifort-debug'
                rm -rf $folder && mkdir $folder && pushd $folder
                FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Debug  -DWITH_TESTING=ON
                make
                ctest --output-on-failure
                '''
            }
        }
        stage('Build and Test Intel-Release') {
            steps {
                sh '''
                module load cmake
                module load intel
                folder='test-ifort'
                rm -rf $folder && mkdir $folder && pushd $folder
                FC=ifort cmake ../  -DCMAKE_BUILD_TYPE=Release  -DWITH_TESTING=ON
                make
                ctest --output-on-failure
                '''
            }
        }        
    }
}
