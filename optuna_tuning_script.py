import optuna
from GNN4 import main
import matplotlib.pyplot as plt
import json

def objective(trial):
    # Suggest hyperparameters
    
    hidden_dim = trial.suggest_categorical("hidden_dim", [32, 64])
    dropout_prob = trial.suggest_float("dropout_prob", 0.2, 0.5)
    learning_rate = trial.suggest_float("learning_rate", 1e-5, 1e-3, log=True)  
    heads = trial.suggest_categorical("heads_user", [2, 4, 6])                 
    num_epochs = trial.suggest_categorical("num_epochs", [50, 100])
    bsu = trial.suggest_categorical("bsu", [16, 32, 64])
    hops = trial.suggest_categorical("hops_user", [1, 2, 3])                   
    optuna_arg = "True"
    # Fixed parameters
    directory = "/Users/harrytriantafyllidis/HARMONIZER/RUN_TESTS/opt_TCGA"
    csv_file = "/Users/harrytriantafyllidis/HARMONIZER/RUN_TESTS/GNN_labels_TCGA.csv"
    tts = 0.7  # Train-test split ratio
    min_obs = 1

    # Run main pipeline
    try:
        misclassification_percentage = main(
            directory=directory,
            csv_file=csv_file,
            num_epochs=num_epochs,
            learning_rate=learning_rate,
            tts=tts,
            min_obs=min_obs,
            bsu=bsu,
            hidden_dim=hidden_dim,
            dropout_prob=dropout_prob,
            heads_user=heads,
            hops_user = hops,
            optuna = optuna_arg
        )
        if misclassification_percentage is None:
            raise ValueError("Returned misclassification is None")
    except Exception as e:
        print("Trial failed:", e)
        return float("inf")

    return float(misclassification_percentage)


def plot_static_performance(study):
    values = [trial.value for trial in study.trials if trial.value is not None and trial.state == optuna.trial.TrialState.COMPLETE]
    plt.figure(figsize=(10, 6))
    plt.plot(values, marker='o', label='Misclassification % per Trial')
    best_value = study.best_value
    plt.axhline(best_value, color='r', linestyle='--', label=f'Best: {best_value:.2f}%')
    plt.title("Misclassification Percentage per Trial")
    plt.xlabel("Trial")
    plt.ylabel("Misclassification %")
    plt.legend()
    plt.grid(True)
    plt.savefig("optuna_trials_performance.png", dpi=600, bbox_inches='tight')
    plt.show()


# Create and run Optuna study
if __name__ == "__main__":

    def early_stop_callback(study, trial):
        if study.best_value is not None and study.best_value < 10.0:
            print("Early stopping: misclassification below 10%")
            return True
        return False

    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=50, callbacks=[early_stop_callback])

    print("Best trial:")
    print(study.best_trial.params)
    
  
    # Save best parameters to JSON
    with open("best_optuna_params.json", "w") as f:
        json.dump(study.best_params, f, indent=4)
    
    print("Best parameters saved to best_optuna_params.json")

    # Static performance plot
    plot_static_performance(study)
