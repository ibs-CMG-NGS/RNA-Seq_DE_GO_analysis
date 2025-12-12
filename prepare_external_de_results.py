#!/usr/bin/env python3
"""
외부에서 수행한 DE 분석 결과를 파이프라인의 downstream 분석에 사용할 수 있도록 준비하는 스크립트

Usage:
    python prepare_external_de_results.py --input your_de_results.csv --comparison H2O2_vs_Control --output-dir output/H2O2_Neuron

Requirements:
    - 입력 파일에 최소한 'symbol', 'log2FoldChange', 'padj' 컬럼이 있어야 함
    - config 파일의 pairwise_comparisons에 해당 비교쌍이 정의되어 있어야 함
"""

import argparse
import pandas as pd
import yaml
from pathlib import Path
import shutil
import sys

def validate_de_results(df, required_columns=None):
    """DE 결과 파일의 필수 컬럼 확인"""
    if required_columns is None:
        required_columns = ['symbol', 'log2FoldChange', 'padj']
    
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        print(f"❌ 오류: 다음 필수 컬럼이 없습니다: {', '.join(missing_cols)}")
        print(f"현재 컬럼: {', '.join(df.columns)}")
        return False
    
    print(f"✓ 필수 컬럼 확인 완료: {', '.join(required_columns)}")
    
    # 추가 권장 컬럼 확인
    recommended_columns = ['baseMean', 'pvalue', 'lfcSE', 'stat']
    available_recommended = [col for col in recommended_columns if col in df.columns]
    if available_recommended:
        print(f"✓ 권장 컬럼 포함: {', '.join(available_recommended)}")
    
    missing_recommended = [col for col in recommended_columns if col not in df.columns]
    if missing_recommended:
        print(f"ℹ️  권장 컬럼 누락 (선택사항): {', '.join(missing_recommended)}")
    
    return True

def prepare_de_file(input_path, comparison, output_dir, config_path=None):
    """
    DE 결과 파일을 파이프라인 형식에 맞게 준비
    
    Args:
        input_path: 외부 DE 결과 파일 경로
        comparison: 비교쌍 이름 (예: "H2O2_vs_Control")
        output_dir: 출력 디렉토리 (예: "output/H2O2_Neuron")
        config_path: config 파일 경로 (검증용, 선택사항)
    """
    print(f"\n{'='*60}")
    print(f"DE 결과 파일 준비 시작")
    print(f"{'='*60}\n")
    
    # 1. 입력 파일 읽기
    print(f"1. 입력 파일 읽기: {input_path}")
    try:
        df = pd.read_csv(input_path)
        print(f"   - 행 수: {len(df):,}")
        print(f"   - 컬럼 수: {len(df.columns)}")
    except Exception as e:
        print(f"❌ 파일 읽기 오류: {e}")
        return False
    
    # 2. 필수 컬럼 검증
    print(f"\n2. 데이터 검증")
    if not validate_de_results(df):
        return False
    
    # 3. config 파일에서 비교쌍 확인 (선택사항)
    if config_path and Path(config_path).exists():
        print(f"\n3. Config 파일 검증: {config_path}")
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            comparisons = config.get('de_analysis', {}).get('pairwise_comparisons', [])
            comparison_strings = [f"{comp[0]}_vs_{comp[1]}" for comp in comparisons]
            
            if comparison in comparison_strings:
                print(f"   ✓ 비교쌍 '{comparison}'이 config에 정의되어 있습니다")
            else:
                print(f"   ⚠️  경고: 비교쌍 '{comparison}'이 config에 없습니다")
                print(f"   정의된 비교쌍: {', '.join(comparison_strings)}")
                response = input("   계속하시겠습니까? (y/n): ")
                if response.lower() != 'y':
                    return False
        except Exception as e:
            print(f"   ⚠️  Config 파일 읽기 오류 (계속 진행): {e}")
    
    # 4. 출력 디렉토리 생성
    output_path = Path(output_dir) / "pairwise" / comparison
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"\n4. 출력 디렉토리 생성: {output_path}")
    
    # 5. 파일 저장
    output_file = output_path / "final_de_results.csv"
    print(f"\n5. 파일 저장: {output_file}")
    
    # 인덱스 열 처리 (파이프라인 형식에 맞춤)
    # 첫 번째 열이 유전자 ID인 경우 인덱스로 설정
    if df.columns[0] in ['gene_id', 'ensembl_id', 'gene', '']:
        df.set_index(df.columns[0], inplace=True)
    
    df.to_csv(output_file)
    print(f"   ✓ 파일 저장 완료")
    
    # 6. 요약 정보
    print(f"\n{'='*60}")
    print(f"✅ 준비 완료!")
    print(f"{'='*60}")
    print(f"\n다음 단계:")
    print(f"1. Enrichment 분석 실행:")
    print(f"   snakemake --cores 4 {output_dir}/pairwise/{comparison}/.enrichment_done.flag")
    print(f"\n2. 전체 downstream 분석 실행:")
    print(f"   snakemake --cores 4 {output_dir}/pairwise/{comparison}/final_go_results.xlsx")
    print(f"\n3. 특정 분석만 실행:")
    print(f"   snakemake --cores 4 --forcerun go_enrichment kegg_enrichment")
    print()
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='외부 DE 결과를 파이프라인 downstream 분석용으로 준비',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
예시:
  # 단일 비교쌍
  python prepare_external_de_results.py \\
      --input my_de_results.csv \\
      --comparison H2O2_vs_Control \\
      --output-dir output/H2O2_Neuron \\
      --config config_H2O2_Neuron.yml

  # Config 없이 실행
  python prepare_external_de_results.py \\
      --input my_de_results.csv \\
      --comparison GABA_vs_Control \\
      --output-dir output/MyExperiment

필수 컬럼:
  - symbol: 유전자 심볼
  - log2FoldChange: log2 fold change
  - padj: Adjusted p-value

권장 컬럼:
  - baseMean: 평균 발현량
  - pvalue: Raw p-value
  - lfcSE: Log fold change standard error
  - stat: Test statistic
        """
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='입력 DE 결과 CSV 파일 경로')
    parser.add_argument('-c', '--comparison', required=True,
                       help='비교쌍 이름 (예: H2O2_vs_Control, GABA_vs_Control)')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='출력 디렉토리 (예: output/H2O2_Neuron)')
    parser.add_argument('--config', default=None,
                       help='Config YAML 파일 경로 (검증용, 선택사항)')
    parser.add_argument('--required-columns', nargs='+',
                       default=['symbol', 'log2FoldChange', 'padj'],
                       help='필수 컬럼 목록 (기본값: symbol log2FoldChange padj)')
    
    args = parser.parse_args()
    
    # 입력 파일 존재 확인
    if not Path(args.input).exists():
        print(f"❌ 오류: 입력 파일을 찾을 수 없습니다: {args.input}")
        sys.exit(1)
    
    # 실행
    success = prepare_de_file(
        input_path=args.input,
        comparison=args.comparison,
        output_dir=args.output_dir,
        config_path=args.config
    )
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
