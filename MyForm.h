#pragma once
#include "Dirikhle.h"
#include <thread>
#include <msclr\marshal_cppstd.h>

namespace DirikhleGUI {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	using namespace System::Threading;
	using namespace System::Diagnostics;
	using namespace System::ComponentModel;

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^ label1;
	protected:
	private: System::Windows::Forms::TextBox^ textBox1;
	private: System::Windows::Forms::Label^ label2;
	private: System::Windows::Forms::TextBox^ textBox2;
	private: System::Windows::Forms::Label^ label3;
	private: System::Windows::Forms::TextBox^ textBox3;
	private: System::Windows::Forms::Label^ label4;
	private: System::Windows::Forms::TextBox^ textBox4;
	private: System::Windows::Forms::PictureBox^ pictureBox1;
	private: System::Windows::Forms::Label^ label5;
	private: System::Windows::Forms::TextBox^ textBox5;
	private: System::Windows::Forms::DataGridView^ dataGridView1;
	private: System::Windows::Forms::Button^ button1;
	private: System::ComponentModel::IContainer^ components;

	private: System::Windows::Forms::Label^ label7;
	private: System::Windows::Forms::Button^ button2;
	private: System::Windows::Forms::Label^ label6;





	private: Thread^ thread2 = nullptr;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::ComponentModel::ComponentResourceManager^ resources = (gcnew System::ComponentModel::ComponentResourceManager(MyForm::typeid));
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->button2 = (gcnew System::Windows::Forms::Button());
			this->label6 = (gcnew System::Windows::Forms::Label());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(12, 168);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(121, 13);
			this->label1->TabIndex = 0;
			this->label1->Text = L"Число разбиений по X";
			this->label1->Click += gcnew System::EventHandler(this, &MyForm::label1_Click);
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(172, 165);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(100, 20);
			this->textBox1->TabIndex = 1;
			this->textBox1->Text = L"10";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(12, 194);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(121, 13);
			this->label2->TabIndex = 0;
			this->label2->Text = L"Число разбиений по Y";
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(172, 191);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(100, 20);
			this->textBox2->TabIndex = 1;
			this->textBox2->Text = L"10";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(12, 246);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(107, 13);
			this->label3->TabIndex = 0;
			this->label3->Text = L"Ограничение шагов";
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(172, 243);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(100, 20);
			this->textBox3->TabIndex = 1;
			this->textBox3->Text = L"1000";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(12, 220);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(94, 13);
			this->label4->TabIndex = 0;
			this->label4->Text = L"Точность метода";
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(172, 217);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(100, 20);
			this->textBox4->TabIndex = 1;
			this->textBox4->Text = L"1e-8";
			// 
			// pictureBox1
			// 
			this->pictureBox1->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox1.Image")));
			this->pictureBox1->Location = System::Drawing::Point(15, 12);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(269, 135);
			this->pictureBox1->TabIndex = 2;
			this->pictureBox1->TabStop = false;
			this->pictureBox1->Click += gcnew System::EventHandler(this, &MyForm::pictureBox1_Click);
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(12, 272);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(64, 13);
			this->label5->TabIndex = 0;
			this->label5->Text = L"Параметр t";
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(172, 269);
			this->textBox5->Name = L"textBox5";
			this->textBox5->ReadOnly = true;
			this->textBox5->Size = System::Drawing::Size(100, 20);
			this->textBox5->TabIndex = 1;
			// 
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Location = System::Drawing::Point(326, 12);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->Size = System::Drawing::Size(1212, 651);
			this->dataGridView1->TabIndex = 3;
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(15, 312);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(118, 23);
			this->button1->TabIndex = 4;
			this->button1->Text = L"Старт";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(12, 453);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(0, 13);
			this->label7->TabIndex = 6;
			// 
			// button2
			// 
			this->button2->Location = System::Drawing::Point(87, 269);
			this->button2->Name = L"button2";
			this->button2->Size = System::Drawing::Size(79, 20);
			this->button2->TabIndex = 7;
			this->button2->Text = L"Вычислить";
			this->button2->UseVisualStyleBackColor = true;
			this->button2->Click += gcnew System::EventHandler(this, &MyForm::button2_Click);
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(15, 420);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(70, 13);
			this->label6->TabIndex = 8;
			this->label6->Text = L"Результаты:";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1540, 675);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->button2);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->dataGridView1);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->textBox5);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
		private: void BackgroundWorker1_DoWork() {
			// Sleep 2 seconds to emulate getting data.
			system("python show_plot.py");
		}


	private: System::Void pictureBox1_Click(System::Object^ sender, System::EventArgs^ e) {
	}


	private: System::Void button1_Click(System::Object^ sender, System::EventArgs^ e) {
		int Nmax = Convert::ToDouble(textBox3->Text); // максимальное число итераций (не менее 1)
		int S = 0; // счетчик итераций
		int S_2 = 0; // счетчик итераций
		double eps = Convert::ToDouble(textBox4->Text); // минимально допустимый прирост
		double eps_max = 0; // текущее значение прироста
		double eps_max_2 = 0; // текущее значение прироста
		double eps_cur = 0; // для подсчета текущего значения прироста
		double error_max = 0; // для подсчета текущего значения прироста
		double accuracy = 0; // точность
		double a2, k2, h2; // ненулевые элементы матрицы (-A)

		//func my_func;
		func2 my_func; // Тестовая задача
		const int n = Convert::ToDouble(textBox1->Text); //размерность сетки
		const int m = Convert::ToDouble(textBox2->Text); //размерность сетки
		writeHeader(n, m);

		std::vector<std::vector<double>> v(n + 1); // искомый вектор v
		std::vector<std::vector<double>> v_2(2 * n + 1); // искомый вектор v с половинным шагом
		std::vector<double> r((n - 1) * (m - 1)); // невязка
		double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
		double t = t_optimal((b - a) / n, (d - c) / m, n, m);
		textBox5->Text = ""+(t);
		//if (textBox5->Text != "1.5278640450004206") {
		//	w = std::stod(msclr::interop::marshal_as<std::string>(textBox5->Text));
		//}
		//double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
		h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
		k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
		a2 = -2 * (h2 + k2);
		bool flag = false;

		int i, j; //индексы
		double v_old; // старое значение преобразуемой компоненты вектора v
		double v_new; // новое значение преобразуемой компоненты вектора v

		solve(v_2, my_func, 2 * n, 2 * m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
		S_2 = S;
		S = 0;
		eps_max_2 = eps_max;
		eps_max = 0;
		eps_cur = 0;
		error_max = 0;
		auto max_diff = solve(v, my_func, n, m, a, b, c, d, Nmax, S, eps, eps_max, error_max);


		for (j = 1; j < m; j++)
			for (i = 1; i < n; i++) {
				if (is_inside(i, j, n, m)) {
					double k_1 = h2 * (v[i + 1][j] * ((i + 1 == n) ? 0 : 1)
						+ v[i - 1][j] * ((i - 1 == 0) ? 0 : 1));
					double k_2 = k2 * (v[i][j + 1] * ((j + 1 == m) ? 0 : 1)
						+ v[i][j - 1] * ((j - 1 == 0) ? 0 : 1));
					double k_3 = v[i][j] * a2;
					v_new = (k_3 + k_1 + k_2);
					double p = a + i * (b - a) / n;
					double e = c + j * (d - c) / m;
					double func_t = my_func(p, e);
					double f_t = -(h2 * (v[i + 1][j] * ((i + 1 == n) ? 1 : 0)
						+ v[i - 1][j] * ((i - 1 == 0) ? 1 : 0)) +
						k2 * (v[i][j + 1] * ((j + 1 == m) ? 1 : 0)
							+ v[i][j - 1] * ((j - 1 == 0) ? 1 : 0))) + func_t;
					r[(j - 1) * (n - 1) + (i - 1)] = v_new - f_t;
				}
			}

		double r_norm = 0;
		for (j = 0; j < r.size(); j++)
			if(fabs(r[j])>r_norm)
				r_norm = fabs(r[j]);
		std::ofstream outfile("test.dat");

		for (j = 1; j < m; j++)
			for (i = 1; i < n; i++) {
				if (abs(v[i][j] - v_2[2*i][2*j]) > accuracy) {
					accuracy = abs(v[i][j] - v_2[2*i][2*j]);
				}
			}

		for (i = 0; i < n + 1; i++)
			for (j = 0; j < m + 1; j++) {
				v_new = v[i][j];
				double p = a + i * (b - a) / n;
				double e = c + j * (d - c) / m;
				outfile << p << "\t" << e << "\t" << v_new << "\n";
			}
		//writeTable(n, m, vec_u, a, b, c, d, 8);
		std::string results = writeFinalTable(n, m, v, a, b, c, d, S, S_2, eps, eps_max, eps_max_2, error_max, r_norm, t, accuracy);

		String^ text = gcnew String(results.c_str());

		label7->Text = text+L"\nmax_diff="+(max_diff);

		Process^ myProcess = gcnew Process();
		myProcess->StartInfo->FileName = "python";
		myProcess->StartInfo->Arguments = "show_plot.py";
		myProcess->Start();


		dataGridView1->Rows->Clear();
		dataGridView1->ColumnCount = n+1;
		dataGridView1->RowHeadersWidth = 100;
		dataGridView1->TopLeftHeaderCell->Value = "       Y\\X       ";
		for (int r = 0; r < m+1; ++r) {

			DataGridViewRow^ row = gcnew DataGridViewRow();
			row->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());
			//row->CreateCells(this->dataGridView1);

			for (int c = 0; c < n+1; ++c) {
				DataGridViewCell^ cel = gcnew DataGridViewTextBoxCell();
				cel->Value = setPresision(v[c][r], 4);
				row->Cells->Add(cel);
				dataGridView1->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
			}

			this->dataGridView1->Rows->Add(row);
		}
	}
	private: System::Void button2_Click(System::Object^ sender, System::EventArgs^ e) {
		const int n = Convert::ToDouble(textBox1->Text); //размерность сетки
		const int m = Convert::ToDouble(textBox2->Text); //размерность сетки
		double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
		double t = t_optimal((b - a) / n, (d - c) / m, n, m);
		textBox5->Text = gcnew String(doubleToString(t).c_str());
	}
private: System::Void label1_Click(System::Object^ sender, System::EventArgs^ e) {
}
};
}
