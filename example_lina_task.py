import numpy as np
from matplotlib import pyplot
from rdkit import Chem
from rdkit.Chem import Draw


class MoleculeSpectra:
    def __init__(self):
        self.molecules = {
            'бифенил': 'C1=CC=C(C=C1)C2=CC=CC=C2',
            'терфенил': 'C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=CC=C3',
            'дифенил': 'C1=CC=C(C=C1)C2=CC=CC=C2'  # дифенил - синоним бифенила
        }

    def get_uv_spectrum(self, molecule_name):
        """Генерация модельного УФ спектра"""
        # Параметры для разных молекул
        uv_params = {
            'бифенил': {'peaks': [(250, 18000), (280, 12000), (320, 8000)], 'width': 15},
            'терфенил': {'peaks': [(260, 22000), (300, 15000), (350, 10000)], 'width': 12},
            'дифенил': {'peaks': [(250, 18000), (280, 12000), (320, 8000)], 'width': 15}
        }

        wavelength = np.linspace(200, 400, 500)
        absorbance = np.zeros_like(wavelength)

        params = uv_params.get(molecule_name, uv_params['бифенил'])

        for peak_wl, intensity in params['peaks']:
            absorbance += intensity * np.exp(-((wavelength - peak_wl) / params['width']) ** 2)

        return wavelength, absorbance

    def get_ir_spectrum(self, molecule_name):
        """Генерация модельного ИК спектра"""
        # Параметры для разных молекул
        ir_params = {
            'бифенил': {
                'peaks': [(700, 0.8), (750, 0.9), (1000, 0.6), (1200, 0.7),
                          (1500, 0.9), (1600, 1.0), (3000, 0.8), (3050, 0.7)]
            },
            'терфенил': {
                'peaks': [(700, 0.7), (800, 0.8), (1000, 0.5), (1200, 0.6),
                          (1500, 0.8), (1600, 0.9), (3000, 0.7), (3050, 0.6)]
            },
            'дифенил': {
                'peaks': [(700, 0.8), (750, 0.9), (1000, 0.6), (1200, 0.7),
                          (1500, 0.9), (1600, 1.0), (3000, 0.8), (3050, 0.7)]
            }
        }

        wavenumber = np.linspace(500, 3500, 1000)
        transmittance = np.ones_like(wavenumber)

        params = ir_params.get(molecule_name, ir_params['бифенил'])

        for peak_wn, intensity in params['peaks']:
            transmittance -= intensity * np.exp(-((wavenumber - peak_wn) / 20) ** 2)

        transmittance = np.clip(transmittance, 0, 1)
        return wavenumber, transmittance

    def draw_molecule(self, molecule_name):
        """Отрисовка молекулы"""
        smiles = self.molecules.get(molecule_name)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                return img
        return None

    def plot_spectra(self, molecule_name):
        """Построение графиков спектров"""
        if molecule_name not in self.molecules:
            print(f"Молекула {molecule_name} не найдена. Доступные молекулы: {list(self.molecules.keys())}")
            return

        fig, ((ax1, ax2), (ax3, ax4)) = pyplot.subplots(2, 2, figsize=(15, 12))

        # УФ спектр
        wavelength, absorbance = self.get_uv_spectrum(molecule_name)
        ax1.plot(wavelength, absorbance, 'b-', linewidth=2)
        ax1.set_xlabel('Длина волны (нм)')
        ax1.set_ylabel('Молярный коэффициент поглощения')
        ax1.set_title(f'УФ спектр {molecule_name}')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(200, 400)

        # ИК спектр
        wavenumber, transmittance = self.get_ir_spectrum(molecule_name)
        ax2.plot(wavenumber, transmittance * 100, 'r-', linewidth=2)
        ax2.set_xlabel('Волновое число (см⁻¹)')
        ax2.set_ylabel('Пропускание (%)')
        ax2.set_title(f'ИК спектр {molecule_name}')
        ax2.grid(True, alpha=0.3)


        ax2.set_xlim(500, 3500)
        ax2.invert_xaxis()  # ИК спектры обычно отображаются справа налево

        # Структура молекулы
        img = self.draw_molecule(molecule_name)
        if img is not None:
            ax3.imshow(img)
            ax3.axis('off')
            ax3.set_title(f'Структура {molecule_name}')

        # Комбинированный график
        ax4.plot(wavelength, absorbance, 'b-', linewidth=2, label='УФ спектр')
        ax4_twin = ax4.twinx()
        ax4_twin.plot(wavenumber, transmittance * 100, 'r-', linewidth=2, label='ИК спектр')
        ax4.set_xlabel('Длина волны (нм) / Волновое число (см⁻¹)')
        ax4.set_ylabel('Молярный коэффициент поглощения', color='b')
        ax4_twin.set_ylabel('Пропускание (%)', color='r')
        ax4.set_title(f'Комбинированные спектры {molecule_name}')
        ax4.grid(True, alpha=0.3)
        ax4.set_xlim(200, 400)

        # Добавляем легенду
        lines1, labels1 = ax4.get_legend_handles_labels()
        lines2, labels2 = ax4_twin.get_legend_handles_labels()
        ax4.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        pyplot.tight_layout()
        pyplot.show()


    def interactive_plot(self):
        """Интерактивный выбор молекулы"""
        print("Доступные молекулы:")
        for i, mol_name in enumerate(self.molecules.keys(), 1):
            print(f"{i}. {mol_name}")

        while True:
            try:
                choice = input("\nВыберите молекулу (номер или название, или 'q' для выхода): ")
                if choice.lower() == 'q':
                    break

                if choice.isdigit():
                    mol_list = list(self.molecules.keys())
                    if 1 <= int(choice) <= len(mol_list):
                        molecule_name = mol_list[int(choice) - 1]
                    else:
                        print("Неверный номер")
                        continue
                else:
                    if choice in self.molecules:
                        molecule_name = choice
                    else:
                        print("Молекула не найдена")
                        continue

                self.plot_spectra(molecule_name)

            except KeyboardInterrupt:
                break
            except Exception as e:
                print(f"Ошибка: {e}")


# Дополнительная функция для сравнения нескольких молекул
def compare_molecules():
    """Сравнение спектров нескольких молекул"""
    spectra = MoleculeSpectra()

    fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(15, 6))

    molecules = ['бифенил', 'терфенил']
    colors = ['blue', 'red', 'green']

    # Сравнение УФ спектров
    for i, mol in enumerate(molecules):
        wavelength, absorbance = spectra.get_uv_spectrum(mol)
        ax1.plot(wavelength, absorbance, color=colors[i], linewidth=2, label=mol)

    ax1.set_xlabel('Длина волны (нм)')
    ax1.set_ylabel('Молярный коэффициент поглощения')
    ax1.set_title('Сравнение УФ спектров')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Сравнение ИК спектров
    for i, mol in enumerate(molecules):
        wavenumber, transmittance = spectra.get_ir_spectrum(mol)
        ax2.plot(wavenumber, transmittance * 100, color=colors[i], linewidth=2, label=mol)

    ax2.set_xlabel('Волновое число (см⁻¹)')
    ax2.set_ylabel('Пропускание (%)')
    ax2.set_title('Сравнение ИК спектров')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.invert_xaxis()

    pyplot.tight_layout()
    pyplot.show()


# Использование
if __name__ == "__main__":
    # Создаем экземпляр класса
    spectra_analyzer = MoleculeSpectra()

    # Способ 1: Построить спектр для конкретной молекулы
    print("Построение спектров для бифенила:")
    spectra_analyzer.plot_spectra('бифенил')

    # Способ 2: Интерактивный режим
    print("\nИнтерактивный режим:")
    spectra_analyzer.interactive_plot()

    # Способ 3: Сравнение молекул
    print("\nСравнение молекул:")
    compare_molecules()