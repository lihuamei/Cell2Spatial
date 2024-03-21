import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras import regularizers
from keras.layers import BatchNormalization
from keras.callbacks import ModelCheckpoint, ReduceLROnPlateau, EarlyStopping

def removeCacheModel(directory, prefix):
    files = os.listdir(directory)
    for file in files:
        if not file.startswith(prefix): continue
        file_path = os.path.join(directory, file)
        os.remove(file_path)

def createNetMmodel(out_layers, input_shape):
    model = Sequential()
    model.add(Dense(64, activation='relu', input_shape=(input_shape, ), kernel_regularizer=regularizers.l2(0.01)))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    model.add(Dense(128, activation='relu', kernel_regularizer=regularizers.l2(0.01)))  # Additional hidden layer
    model.add(BatchNormalization())
    model.add(Dropout(0.2))
    model.add(Dense(128, activation='relu', kernel_regularizer=regularizers.l2(0.01)))
    model.add(BatchNormalization())
    model.add(Dense(out_layers, activation='softmax', kernel_regularizer=regularizers.l2(0.01)))
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

def trainModel(X, y_encoded, epochs = 100):
    best_model, best_accuracy = None, 0
    scaler = StandardScaler()
    kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    for fold, (train_index, test_index) in enumerate(kfold.split(X, y_encoded)):
        X_train, X_test = X[train_index], X[test_index]
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        y_train, y_test = y_encoded[train_index], y_encoded[test_index]

        model = createNetMmodel(out_layers = np.unique(y_encoded).shape[0], input_shape = X.shape[1])
        checkpoint = ModelCheckpoint(f'best_model_fold_{fold}.h5', monitor='val_accuracy', save_best_only=True, mode='max',verbose=1)
        reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=5, min_lr=1e-6)
        early_stop = EarlyStopping(monitor='val_accuracy', patience=10, restore_best_weights=True)
        model.fit(X_train, y_train, epochs=epochs, batch_size=32, validation_data=(X_test, y_test), callbacks=[checkpoint, reduce_lr, early_stop])

        scores = model.evaluate(X_test, y_test)
        accuracy = scores[1]
        print(f'Fold {fold + 1} - Accuracy: {accuracy*100:.2f}%')

        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_model = model
    print(f'Best Model Accuracy: {best_accuracy*100:.2f}%')
    return best_model

def runNetModel(sc_data, st_data, st_label, epochs = 100):
    sc_data, st_data, st_label, epochs = np.array(sc_data), np.array(st_data), np.array(st_label), int(epochs)
    model = trainModel(st_data, st_label, epochs)
    removeCacheModel('./', 'best_model_fold_')
    predictions = model.predict(sc_data)
    return pd.DataFrame(predictions)
