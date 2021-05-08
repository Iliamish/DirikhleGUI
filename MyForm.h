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

	private: System::Windows::Forms::Button^ button1;
	private: System::ComponentModel::IContainer^ components;
	private: System::Windows::Forms::Label^ label6;
	private: System::Windows::Forms::Label^ label7;

	private: System::Windows::Forms::ListBox^ listBox1;
	private: System::Windows::Forms::RadioButton^ radioButton1;
	private: System::Windows::Forms::RadioButton^ radioButton2;
	private: System::Windows::Forms::DataGridView^ dataGridView1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ NumberOfColumn;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ X;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ V;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ UandDouble;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ SubVandSome;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Column1;
	private: System::Windows::Forms::DataGridViewTextBoxColumn^ Column2;











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
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->listBox1 = (gcnew System::Windows::Forms::ListBox());
			this->radioButton1 = (gcnew System::Windows::Forms::RadioButton());
			this->radioButton2 = (gcnew System::Windows::Forms::RadioButton());
			this->dataGridView1 = (gcnew System::Windows::Forms::DataGridView());
			this->NumberOfColumn = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->X = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->V = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->UandDouble = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->SubVandSome = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column1 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
			this->Column2 = (gcnew System::Windows::Forms::DataGridViewTextBoxColumn());
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
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(12, 272);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(69, 13);
			this->label5->TabIndex = 0;
			this->label5->Text = L"Параметр w";
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(150, 269);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(100, 20);
			this->textBox5->TabIndex = 1;
			this->textBox5->Text = L"1.5278640450004206";
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
			// dataGridView1
			// 
			this->dataGridView1->ColumnHeadersHeightSizeMode = System::Windows::Forms::DataGridViewColumnHeadersHeightSizeMode::AutoSize;
			this->dataGridView1->Columns->AddRange(gcnew cli::array< System::Windows::Forms::DataGridViewColumn^  >(7) {
				this->NumberOfColumn,
					this->X, this->V, this->UandDouble, this->SubVandSome, this->Column1, this->Column2
			});
			this->dataGridView1->Location = System::Drawing::Point(326, 24);
			this->dataGridView1->Name = L"dataGridView1";
			this->dataGridView1->RowHeadersVisible = false;
			this->dataGridView1->Size = System::Drawing::Size(1132, 612);
			this->dataGridView1->TabIndex = 11;
			this->dataGridView1->CellContentClick += gcnew System::Windows::Forms::DataGridViewCellEventHandler(this, &MyForm::dataGridView1_CellContentClick);
			// 
			// NumberOfColumn
			// 
			this->NumberOfColumn->HeaderText = L"i";
			this->NumberOfColumn->Name = L"NumberOfColumn";
			this->NumberOfColumn->ReadOnly = true;
			// 
			// X
			// 
			this->X->HeaderText = L"Xi-1";
			this->X->Name = L"X";
			this->X->ReadOnly = true;
			this->X->Width = 50;
			// 
			// V
			// 
			this->V->HeaderText = L"Xi";
			this->V->Name = L"V";
			this->V->ReadOnly = true;
			// 
			// UandDouble
			// 
			this->UandDouble->HeaderText = L"ai";
			this->UandDouble->Name = L"UandDouble";
			this->UandDouble->ReadOnly = true;
			// 
			// SubVandSome
			// 
			this->SubVandSome->HeaderText = L"bi";
			this->SubVandSome->Name = L"SubVandSome";
			this->SubVandSome->ReadOnly = true;
			this->SubVandSome->Width = 124;
			// 
			// Column1
			// 
			this->Column1->HeaderText = L"ci";
			this->Column1->Name = L"Column1";
			// 
			// Column2
			// 
			this->Column2->HeaderText = L"di";
			this->Column2->Name = L"Column2";
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(1540, 648);
			this->Controls->Add(this->dataGridView1);
			this->Controls->Add(this->radioButton2);
			this->Controls->Add(this->radioButton1);
			this->Controls->Add(this->listBox1);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->button1);
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
		double accuracy = 1000; // точность
		double a2, k2, h2; // ненулевые элементы матрицы (-A)


		const int n = Convert::ToDouble(textBox1->Text); //размерность сетки
		const int m = Convert::ToDouble(textBox2->Text); //размерность сетки
		writeHeader(n, m);
		std::vector<std::vector<double> >  v(n+1, std::vector<double>(m+1));
		//std::vector<std::vector<double>> v(n + 1); // искомый вектор v
		std::vector<std::vector<double> >  v_2(2*n + 1, std::vector<double>(2*m + 1)); // искомый вектор v с половинным шагом
		std::vector<double> r((n - 1) * (m - 1)); // невязка
		double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
		double w = 1.5278640450004206;
		if (textBox5->Text != "1.5278640450004206") {
			w = std::stod(msclr::interop::marshal_as<std::string>(textBox5->Text));
		}
		func my_func;
		//double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
		h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
		k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
		a2 = 2 * (h2 + k2);
		bool flag = false;

		int i, j; //индексы
		double v_old; // старое значение преобразуемой компоненты вектора v
		double v_new; // новое значение преобразуемой компоненты вектора v
		a = -1;
		b = 1;
		TResults result= tfunc1(a, b, n, testFunc);
	
		for (size_t i = 0; i < result.b.size(); i++) {
	
			dataGridView1->Rows->Add();
			dataGridView1->Rows[i]->Cells[0]->Value = i+1;
			dataGridView1->Rows[i]->Cells[1]->Value = equalsZero(result.res_vec[i].first);
			dataGridView1->Rows[i]->Cells[2]->Value = equalsZero(result.res_vec[i + 1].first);
			dataGridView1->Rows[i]->Cells[3]->Value = equalsZero(result.a[i]);
			dataGridView1->Rows[i]->Cells[4]->Value = equalsZero(result.b[i]);
			dataGridView1->Rows[i]->Cells[5]->Value = equalsZero(result.c[i]);
			dataGridView1->Rows[i]->Cells[6]->Value = equalsZero (result.d[i]);
		}
		
	}
private: System::Void checkBox1_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void listBox1_SelectedIndexChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void radioButton2_CheckedChanged(System::Object^ sender, System::EventArgs^ e) {
}
private: System::Void dataGridView1_CellContentClick(System::Object^ sender, System::Windows::Forms::DataGridViewCellEventArgs^ e) {
}
};
}
