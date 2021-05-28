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


	private: System::Windows::Forms::DataGridView^ dataGridView1;
	private: System::Windows::Forms::Button^ button1;
	private: System::ComponentModel::IContainer^ components;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::Label^ label7;

	private: System::Windows::Forms::ListBox^ listBox1;
	private: System::Windows::Forms::RadioButton^ radioButton1;
	private: System::Windows::Forms::RadioButton^ radioButton2;
	private: System::Windows::Forms::TabControl^ tabControl1;
	private: System::Windows::Forms::TabPage^ tabPage1;
	private: System::Windows::Forms::TabPage^ tabPage2;
	private: System::Windows::Forms::DataGridView^ dataGridView2;
	private: System::Windows::Forms::TabPage^ tabPage3;
	private: System::Windows::Forms::DataGridView^ dataGridView3;

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
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->listBox1 = (gcnew System::Windows::Forms::ListBox());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->tabControl1 = (gcnew System::Windows::Forms::TabControl());
			this->tabPage1 = (gcnew System::Windows::Forms::TabPage());
			this->tabPage2 = (gcnew System::Windows::Forms::TabPage());
			this->dataGridView2 = (gcnew System::Windows::Forms::DataGridView());
			this->tabPage3 = (gcnew System::Windows::Forms::TabPage());
			this->dataGridView3 = (gcnew System::Windows::Forms::DataGridView());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView1))->BeginInit();
			this->tabControl1->SuspendLayout();
			this->tabPage1->SuspendLayout();
			this->tabPage2->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->BeginInit();
			this->tabPage3->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView3))->BeginInit();
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
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(150, 165);
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
			this->textBox2->Location = System::Drawing::Point(150, 191);
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
			this->textBox3->Location = System::Drawing::Point(150, 243);
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
			this->textBox4->Location = System::Drawing::Point(150, 217);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(100, 20);
			this->textBox4->TabIndex = 1;
			this->textBox4->Text = L"1e-8";
			// 
			// pictureBox1
			// 
			this->pictureBox1->Image = (cli::safe_cast<System::Drawing::Image^>(resources->GetObject(L"pictureBox1.Image")));
			this->pictureBox1->Location = System::Drawing::Point(15, 24);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(269, 135);
			this->pictureBox1->TabIndex = 2;
			this->pictureBox1->TabStop = false;
			this->pictureBox1->Click += gcnew System::EventHandler(this, &MyForm::pictureBox1_Click);
			// 
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Location = System::Drawing::Point(-4, -4);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->Size = System::Drawing::Size(1212, 601);
			this->dataGridView1->TabIndex = 3;
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(12, 324);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(118, 23);
			this->button1->TabIndex = 4;
			this->button1->Text = L"Вычислить";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(15, 361);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(70, 13);
			this->label6->TabIndex = 5;
			this->label6->Text = L"Результаты:";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(15, 386);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(0, 13);
			this->label7->TabIndex = 6;
			// 
			// listBox1
			// 
			this->listBox1->FormattingEnabled = true;
			this->listBox1->Items->AddRange(gcnew cli::array< System::Object^  >(5) {
				L"Метод Верхней Релаксации", L"Сопряженные градиенты (RECT)",
					L"Сопряженные градиенты (CUT)", L"Метод Простых Итераций (RECT)", L"Метод Простых Итераций (CUT)"
			});
			this->listBox1->Location = System::Drawing::Point(136, 324);
			this->listBox1->Name = L"listBox1";
			this->listBox1->Size = System::Drawing::Size(184, 30);
			this->listBox1->TabIndex = 8;
			this->listBox1->SelectedIndexChanged += gcnew System::EventHandler(this, &MyForm::listBox1_SelectedIndexChanged);
			// 
			// radioButton1
			// 
			this->radioButton1->AutoSize = true;
			this->radioButton1->Location = System::Drawing::Point(15, 301);
			this->radioButton1->Name = L"radioButton1";
			this->radioButton1->Size = System::Drawing::Size(73, 17);
			this->radioButton1->TabIndex = 9;
			this->radioButton1->TabStop = true;
			this->radioButton1->Text = L"Тестовая";
			this->radioButton1->UseVisualStyleBackColor = true;
			// 
			// radioButton2
			// 
			this->radioButton2->AutoSize = true;
			this->radioButton2->Location = System::Drawing::Point(128, 301);
			this->radioButton2->Name = L"radioButton2";
			this->radioButton2->Size = System::Drawing::Size(75, 17);
			this->radioButton2->TabIndex = 10;
			this->radioButton2->TabStop = true;
			this->radioButton2->Text = L"Основная";
			this->radioButton2->UseVisualStyleBackColor = true;
			this->radioButton2->CheckedChanged += gcnew System::EventHandler(this, &MyForm::radioButton2_CheckedChanged);
			// 
			// tabControl1
			// 
			this->tabControl1->Controls->Add(this->tabPage1);
			this->tabControl1->Controls->Add(this->tabPage2);
			this->tabControl1->Controls->Add(this->tabPage3);
			this->tabControl1->Location = System::Drawing::Point(326, 13);
			this->tabControl1->Name = L"tabControl1";
			this->tabControl1->SelectedIndex = 0;
			this->tabControl1->Size = System::Drawing::Size(1202, 623);
			this->tabControl1->TabIndex = 11;
			// 
			// tabPage1
			// 
			this->tabPage1->Controls->Add(this->dataGridView1);
			this->tabPage1->Location = System::Drawing::Point(4, 22);
			this->tabPage1->Name = L"tabPage1";
			this->tabPage1->Padding = System::Windows::Forms::Padding(3);
			this->tabPage1->Size = System::Drawing::Size(1194, 597);
			this->tabPage1->TabIndex = 0;
			this->tabPage1->Text = L"V(x,y)";
			this->tabPage1->UseVisualStyleBackColor = true;
			// 
			// tabPage2
			// 
			this->tabPage2->Controls->Add(this->dataGridView2);
			this->tabPage2->Location = System::Drawing::Point(4, 22);
			this->tabPage2->Name = L"tabPage2";
			this->tabPage2->Padding = System::Windows::Forms::Padding(3);
			this->tabPage2->Size = System::Drawing::Size(1194, 597);
			this->tabPage2->TabIndex = 1;
			this->tabPage2->Text = L"U(x,y)";
			this->tabPage2->UseVisualStyleBackColor = true;
			// 
			// dataGridView2
			// 
			this->dataGridView2->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView2->Location = System::Drawing::Point(-4, -4);
			this->dataGridView2->Name = L"dataGridView2";
			this->dataGridView2->Size = System::Drawing::Size(1212, 601);
			this->dataGridView2->TabIndex = 4;
			// 
			// tabPage3
			// 
			this->tabPage3->Controls->Add(this->dataGridView3);
			this->tabPage3->Location = System::Drawing::Point(4, 22);
			this->tabPage3->Name = L"tabPage3";
			this->tabPage3->Padding = System::Windows::Forms::Padding(3);
			this->tabPage3->Size = System::Drawing::Size(1194, 597);
			this->tabPage3->TabIndex = 2;
			this->tabPage3->Text = L"U(x,y) - V(x,y)";
			this->tabPage3->UseVisualStyleBackColor = true;
			// 
			// dataGridView3
			// 
			this->dataGridView3->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView3->Location = System::Drawing::Point(-4, -4);
			this->dataGridView3->Name = L"dataGridView3";
			this->dataGridView3->Size = System::Drawing::Size(1212, 601);
			this->dataGridView3->TabIndex = 5;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1540, 648);
			this->Controls->Add(this->tabControl1);
			this->Controls->Add(this->radioButton2);
			this->Controls->Add(this->radioButton1);
			this->Controls->Add(this->listBox1);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->label4);
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
			this->tabControl1->ResumeLayout(false);
			this->tabPage1->ResumeLayout(false);
			this->tabPage2->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView2))->EndInit();
			this->tabPage3->ResumeLayout(false);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->dataGridView3))->EndInit();
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
		double accuracy = 1000; // точность
		double a2, k2, h2; // ненулевые элементы матрицы (-A)


		const int n = Convert::ToDouble(textBox1->Text); //размерность сетки
		const int m = Convert::ToDouble(textBox2->Text); //размерность сетки
		writeHeader(n, m);
		std::vector<std::vector<double> >  v(n+1, std::vector<double>(m+1));
		std::vector<std::vector<double> >  u(n + 1, std::vector<double>(m + 1));
		//std::vector<std::vector<double>> v(n + 1); // искомый вектор v
		std::vector<std::vector<double> >  v_2(2*n + 1, std::vector<double>(2*m + 1)); // искомый вектор v с половинным шаго
		double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
		double w = 1.5278640450004206;
	
		func my_func;
		//double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
		h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
		k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
		a2 = 2 * (h2 + k2);
		bool flag = false;
		results r;
		int i, j; //индексы
		double v_old; // старое значение преобразуемой компоненты вектора v
		double v_new; // новое значение преобразуемой компоненты вектора v
		if (radioButton1->Checked) {
			if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (RECT)") {
				fillU(u, a, b, c, d, n, m);
			}  else {
				if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (CUT)")
				{
					fill_u_test(u, a, b, c, d, n, m);
				}
		} 
		}
		else {
			if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (RECT)") {
				r =solve(v, funcDef, n, m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
			}
		}

	
		S_2 = 0;
		eps_max_2 = eps_max;
		eps_max = 0;
		eps_cur = 0;
		error_max = 0;


		if (radioButton1->Checked) {
			if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (RECT)") {
				r = solveTest(v, funcTest, n,  m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
			}
			else {
				if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (CUT)")
				{
					r = solveCut(v, funcTest, n, m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
				}
			}
		}
		else {
			 if (listBox1->SelectedItem->ToString() == "Сопряженные градиенты (RECT)") {
				solve(v_2, funcDef,  2*n, 2* m, a, b, c, d, Nmax, S_2, eps, eps_max, error_max);
			}
		}


		
	


		std::ofstream outfile("test.dat");
		accuracy = 0;
		if (radioButton2->Checked) {
			for (j = 1; j < m; j++)
				for (i = 1; i < n; i++) {
					if (abs(v[i][j] - v_2[2 * i][2 * j]) > accuracy) {
						accuracy = abs(v[i][j] - v_2[2*i][2*j]);
					}
				}
		}
		else {
			for (j = 1; j < m; j++)
				for (i = 1; i < n; i++) {
					if (abs(v[i][j] - u[i][j]) > accuracy) {
						accuracy = abs(v[i][j] - u[i][j]);
					}
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
		std::string results = writeFinalTable(n, m, v, a, b, c, d, S, S_2, eps, eps_max, eps_max_2, error_max, 0, w, accuracy,r);

		String^ text = gcnew String(results.c_str());

		label7->Text = text;

		
		Process^ myProcess = gcnew Process();
		myProcess->StartInfo->FileName = "python";
		myProcess->StartInfo->Arguments = "show_plot.py";
		myProcess->Start();

		

		dataGridView1->Rows->Clear();
		dataGridView1->ColumnCount = n+1;
		dataGridView1->RowHeadersWidth = 100;
		dataGridView1->TopLeftHeaderCell->Value = "       Y\\X       ";

		dataGridView2->Rows->Clear();
		dataGridView2->ColumnCount = n + 1;
		dataGridView2->RowHeadersWidth = 100;
		dataGridView2->TopLeftHeaderCell->Value = "       Y\\X       ";

		dataGridView3->Rows->Clear();
		dataGridView3->ColumnCount = n + 1;
		dataGridView3->RowHeadersWidth = 100;
		dataGridView3->TopLeftHeaderCell->Value = "       Y\\X       ";


		if (radioButton1->Checked) {
			for (int r = m; r >= 0; r--) {
				DataGridViewRow^ row1 = gcnew DataGridViewRow();
				row1->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());

				DataGridViewRow^ row2 = gcnew DataGridViewRow();
				row2->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());

				DataGridViewRow^ row3 = gcnew DataGridViewRow();
				row3->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());

				for (int c = 0; c < n + 1; c++) {
					DataGridViewCell^ cel1 = gcnew DataGridViewTextBoxCell();
					DataGridViewCell^ cel2 = gcnew DataGridViewTextBoxCell();
					DataGridViewCell^ cel3 = gcnew DataGridViewTextBoxCell();
					cel1->Value = setPresision(v[c][r], 3);
					row1->Cells->Add(cel1);

					cel2->Value = setPresision(u[c][r], 3);
					row2->Cells->Add(cel2);

					cel3->Value = setPresision(u[c][r] - v[c][r], 3);
					row3->Cells->Add(cel3);

					dataGridView1->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
					dataGridView2->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
					dataGridView3->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
				}

				this->dataGridView1->Rows->Add(row1);
				this->dataGridView2->Rows->Add(row2);
				this->dataGridView3->Rows->Add(row3);
			}
		}
		else {
			dataGridView2->ColumnCount = 2*n + 1;
			for (int r = m; r >= 0; r--) {
				tabControl1->TabPages[1]->Text = "V2(x,y)";
				tabControl1->TabPages[2]->Text = "V2(x,y) - V(x,y)";
				DataGridViewRow^ row1 = gcnew DataGridViewRow();
				row1->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());

				DataGridViewRow^ row3 = gcnew DataGridViewRow();
				row3->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / m).c_str());

				for (int c = 0; c < n + 1; c++) {
					DataGridViewCell^ cel1 = gcnew DataGridViewTextBoxCell();
					DataGridViewCell^ cel3 = gcnew DataGridViewTextBoxCell();
					cel1->Value = setPresision(v[c][r], 3);
					row1->Cells->Add(cel1);


					cel3->Value = setPresision(v_2[2*c][2*r] - v[c][r], 3);
					row3->Cells->Add(cel3);

					dataGridView1->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
					dataGridView3->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / n).c_str());
				}

				this->dataGridView1->Rows->Add(row1);
				this->dataGridView3->Rows->Add(row3);
			}

			for (int r = 2 * m; r >= 0; r--) {
				DataGridViewRow^ row2 = gcnew DataGridViewRow();
				row2->HeaderCell->Value = gcnew String(std::to_string(c + r * (d - c) / (2*m)).c_str());
				for (int c = 0; c < 2*n + 1; c++) {
					DataGridViewCell^ cel2 = gcnew DataGridViewTextBoxCell();
					cel2->Value = setPresision(v_2[c][r], 3);
					row2->Cells->Add(cel2);
					dataGridView2->Columns[c]->Name = gcnew String(std::to_string(a + c * (b - a) / (2*n)).c_str());
				}
				this->dataGridView2->Rows->Add(row2);
			}
			
		}
	

	}
private: System::Void checkBox1_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void listBox1_SelectedIndexChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void radioButton2_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
};
}
