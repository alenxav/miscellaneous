import pip

def install(package):
    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])

install('cython')

install('pip')

install('pandas')
    
install('arrow')

install('numpy')

install('scikit-learn')

install('tensorflow')

install('keras')

install('xgboost')
